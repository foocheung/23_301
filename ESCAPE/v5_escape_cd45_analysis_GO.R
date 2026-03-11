#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(limma)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(glue)
  library(fs)
  library(tidyr)
})

set.seed(42)

logmsg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  cat(msg, "\n")
  msg
}

stage <- function(n, title) {
  banner <- paste0(
    "\n", strrep("=", 70), "\n",
    sprintf(" STAGE %d — %s", n, title), "\n",
    strrep("=", 70), "\n"
  )
  cat(banner)
  invisible(NULL)
}

sanitize_filename <- function(x) gsub("[/\\\\:*?\"<>|]", "_", x)

# =============================================================================
# STAGE 1: Configuration
# =============================================================================
stage(1, "Configuration")

SEURAT_RDS <- "/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/object_final_annotated_panGI_withDay_ESCAPE_Reactome_ALLCELLS_ARRAY.rds"
ESCAPE_REACTOME_RDS <- "/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/ESCAPE_Reactome_UCell_scores_ALLCELLS_ARRAY.rds"
ESCAPE_GOBP_RDS <- "/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/v2_ESCAPE_GO_UCell_scores_ALLCELLS_ARRAY.rds"

SUBJECT_COL <- "Subject_fixed"
TIME_COL    <- "Treatment"   # will be normalized to exactly Untreated / Treated
CD45_COL    <- "CD45pos"

MIN_PATHWAYSUM <- 0.1
MIN_SAMPLES    <- 4L
FDR_CUTOFF     <- 0.05
LOGFC_CUTOFF   <- 0.25

ALLOW_UNPAIRED_FALLBACK <- TRUE
REQUIRE_GOBP <- TRUE   # NEW: fail if GO:BP scores file is missing

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_ROOT  <- paste0("v2_ESCAPE_OUTPUT_FULL_", timestamp)
dir_create(OUT_ROOT)
logmsg("Output root: %s", OUT_ROOT)

DIR_DEG  <- file.path(OUT_ROOT, "deg_results")
DIR_PLOT <- file.path(OUT_ROOT, "plots")
DIR_QC   <- file.path(OUT_ROOT, "qc")
dir_create(c(DIR_DEG, DIR_PLOT, DIR_QC))

# =============================================================================
# STAGE 2: Load data
# =============================================================================
stage(2, "Load data")

logmsg("Loading Seurat object: %s", SEURAT_RDS)
obj <- readRDS(SEURAT_RDS)
logmsg("  Seurat cells: %d", ncol(obj))

logmsg("Loading ESCAPE Reactome scores: %s", ESCAPE_REACTOME_RDS)
escape_reactome <- readRDS(ESCAPE_REACTOME_RDS)
logmsg("  Reactome pathways: %d", nrow(escape_reactome))

escape_gobp <- NULL
if (file.exists(ESCAPE_GOBP_RDS)) {
  logmsg("Loading ESCAPE GO:BP scores: %s", ESCAPE_GOBP_RDS)
  escape_gobp <- readRDS(ESCAPE_GOBP_RDS)
  logmsg("  GO:BP pathways: %d", nrow(escape_gobp))
} else {
  msg <- sprintf("GO:BP file not found at: %s", ESCAPE_GOBP_RDS)
  if (REQUIRE_GOBP) stop(msg)
  logmsg("WARNING: %s", msg)
  logmsg("  GO:BP will be skipped. (Set REQUIRE_GOBP <- TRUE to hard-fail.)")
}

meta_all <- obj@meta.data
rownames(meta_all) <- colnames(obj)

required_cols <- c(SUBJECT_COL, TIME_COL, CD45_COL)
missing <- setdiff(required_cols, colnames(meta_all))
if (length(missing) > 0) stop("Missing metadata columns: ", paste(missing, collapse = ", "))

# ---- NEW: normalize Treatment labels to exactly Untreated / Treated ----
meta_all[[TIME_COL]] <- as.character(meta_all[[TIME_COL]])
meta_all[[TIME_COL]] <- dplyr::case_when(
  meta_all[[TIME_COL]] %in% c("Untreated","Control","NoTx","No_Tx","Baseline","Vehicle","Unstim","Unstimulated","PMA-") ~ "Untreated",
  meta_all[[TIME_COL]] %in% c("Treated","Tx","Treatment","Stim","Stimulated","PMA","PMA+","PMA_Treated") ~ "Treated",
  TRUE ~ meta_all[[TIME_COL]]
)

trt_levels <- sort(unique(meta_all[[TIME_COL]]))
logmsg("Treatment levels after normalization: %s", paste(trt_levels, collapse=", "))
if (!all(c("Untreated","Treated") %in% trt_levels)) {
  stop("Treatment levels must include Untreated and Treated after normalization. Found: ",
       paste(trt_levels, collapse=", "))
}

# =============================================================================
# STAGE 3: Helper functions
# =============================================================================
stage(3, "Helper functions")

pseudobulk_mean <- function(scores, groups) {
  idx <- split(seq_len(ncol(scores)), groups)
  mats <- lapply(idx, function(cols) {
    if (length(cols) == 1) {
      scores[, cols, drop = FALSE]
    } else {
      Matrix::rowMeans(scores[, cols, drop = FALSE])
    }
  })
  m <- do.call(cbind, mats)
  colnames(m) <- names(idx)
  as.matrix(m)
}

add_cd45_status <- function(meta_df) {
  x <- meta_df[[CD45_COL]]
  
  if (is.logical(x)) {
    cd45_log <- x
  } else if (is.numeric(x)) {
    cd45_log <- ifelse(is.na(x), NA, x != 0)
  } else {
    xs <- tolower(trimws(as.character(x)))
    cd45_log <- dplyr::case_when(
      xs %in% c("true","t","1","yes","y","pos","positive","cd45+","cd45pos","cd45_pos","cd45plus","+") ~ TRUE,
      xs %in% c("false","f","0","no","n","neg","negative","cd45-","cd45neg","cd45_neg","cd45minus","-") ~ FALSE,
      TRUE ~ NA
    )
  }
  
  meta_df$immune_status <- dplyr::case_when(
    cd45_log %in% TRUE  ~ "Immune",
    cd45_log %in% FALSE ~ "NonImmune",
    TRUE ~ "Other"
  )
  
  meta_df$CD45pos_coerced <- cd45_log
  meta_df
}

plot_volcano <- function(df, title, out_png) {
  # expects df has logFC + adj.P.Val
  df <- df %>%
    mutate(sig = adj.P.Val <= FDR_CUTOFF & abs(logFC) >= LOGFC_CUTOFF)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(alpha = sig, color = sig)) +
    scale_color_manual(values = c("TRUE"="red", "FALSE"="gray50"), guide="none") +
    scale_alpha_manual(values = c("TRUE"=1, "FALSE"=0.4), guide="none") +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = 2, color = "blue") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = 2, color = "blue") +
    labs(title = title, x = "logFC", y = "-log10(FDR)") +
    theme_minimal(base_size = 11)
  
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

# ---- NEW: small QC writer for pseudobulk build ----
write_build_qc <- function(pathway_type, tag, st, pbk) {
  qc_samples <- st %>%
    mutate(
      n_pathways = nrow(pbk),
      n_samples  = ncol(pbk)
    )
  
  out_qc_samples <- file.path(DIR_QC, glue("{pathway_type}_QC_samples_{sanitize_filename(tag)}.csv"))
  fwrite(qc_samples, out_qc_samples)
  
  out_qc_dims <- file.path(DIR_QC, glue("{pathway_type}_QC_dims_{sanitize_filename(tag)}.txt"))
  txt <- c(
    paste0("tag: ", tag),
    paste0("pathway_type: ", pathway_type),
    paste0("n_pathways: ", nrow(pbk)),
    paste0("n_samples: ", ncol(pbk)),
    paste0("subjects: ", length(unique(st$subject))),
    paste0("immune_status_levels: ", paste(levels(factor(st$immune_status)), collapse=", ")),
    paste0("treatment_levels: ", paste(levels(factor(st$treatment)), collapse=", "))
  )
  writeLines(txt, out_qc_dims)
}

build_pseudobulk <- function(escape_mat, meta_df) {
  cat("\n=== BUILD_PSEUDOBULK ===\n")
  
  meta_df <- add_cd45_status(meta_df)
  meta_df <- meta_df %>% filter(immune_status %in% c("Immune","NonImmune"))
  if (nrow(meta_df) == 0) return(list(error="no_cells_after_cd45_filtering"))
  
  cells_use <- intersect(rownames(meta_df), colnames(escape_mat))
  if (length(cells_use) == 0) return(list(error="no_cellname_overlap"))
  
  meta_df <- meta_df[cells_use, , drop=FALSE]
  
  # pseudobulk key: subject x immune_status x treatment
  meta_df$pb_key <- paste(meta_df[[SUBJECT_COL]], meta_df$immune_status, meta_df[[TIME_COL]], sep="__")
  
  pb <- pseudobulk_mean(escape_mat[, rownames(meta_df), drop=FALSE], meta_df$pb_key)
  
  st <- meta_df %>%
    select(
      pb_key,
      subject = all_of(SUBJECT_COL),
      immune_status,
      treatment = all_of(TIME_COL)
    ) %>%
    distinct()
  
  st <- st[match(colnames(pb), st$pb_key), , drop=FALSE]
  rownames(st) <- st$pb_key
  
  if (ncol(pb) < MIN_SAMPLES) return(list(error=sprintf("few_samples(%d)", ncol(pb))))
  
  pathway_sums <- Matrix::rowSums(pb)
  keep <- pathway_sums >= MIN_PATHWAYSUM
  pbk <- pb[keep, , drop=FALSE]
  
  # Remove pathways with any NA/Inf values
  has_bad <- apply(pbk, 1, function(x) any(is.na(x) | is.infinite(x)))
  pbk <- pbk[!has_bad, , drop=FALSE]
  
  if (nrow(pbk) < 5) return(list(error=sprintf("few_pathways(%d)", nrow(pbk))))
  
  cat("Samples:", ncol(pbk), "| Pathways:", nrow(pbk), "\n")
  list(pbk=pbk, st=st)
}

run_limma_cd45 <- function(pbk, st, allow_fallback=TRUE) {
  immune_status <- factor(st$immune_status, levels = c("NonImmune","Immune"))
  subject <- factor(st$subject)
  
  # prefer paired subjects (have both immune statuses)
  paired_subjects <- st %>%
    distinct(subject, immune_status) %>%
    count(subject) %>%
    filter(n >= 2) %>%
    pull(subject)
  
  st_p <- st %>% filter(subject %in% paired_subjects)
  
  if (nrow(st_p) >= MIN_SAMPLES && length(unique(st_p$subject)) >= 2) {
    pbk_p <- pbk[, rownames(st_p), drop=FALSE]
    immune_status_p <- factor(st_p$immune_status, levels=c("NonImmune","Immune"))
    subject_p <- factor(st_p$subject)
    
    design <- model.matrix(~ subject_p + immune_status_p)
    fit <- lmFit(pbk_p, design)
    fit <- tryCatch(eBayes(fit, trend=TRUE, robust=TRUE), error=function(e) NULL)
    
    if (!is.null(fit)) {
      coef_name <- grep("immune_status", colnames(fit$coefficients), value=TRUE)[1]
      tt <- topTable(fit, coef=coef_name, number=Inf, sort.by="none")
      res <- data.frame(
        pathway = rownames(tt),
        logFC = tt$logFC,
        AveExpr = tt$AveExpr,
        t = tt$t,
        P.Value = tt$P.Value,
        adj.P.Val = tt$adj.P.Val,
        B = tt$B,
        stringsAsFactors = FALSE
      )
      return(list(result=res, mode="paired"))
    }
  }
  
  if (!allow_fallback) return(list(error="no_paired_subjects"))
  
  design2 <- model.matrix(~ immune_status)
  corfit <- duplicateCorrelation(pbk, design2, block=subject)
  fit2 <- lmFit(pbk, design2, block=subject, correlation=corfit$consensus)
  fit2 <- tryCatch(eBayes(fit2, trend=TRUE, robust=TRUE), error=function(e) NULL)
  if (is.null(fit2)) return(list(error="eBayes_failed"))
  
  tt2 <- topTable(fit2, coef="immune_statusImmune", number=Inf, sort.by="none")
  res2 <- data.frame(
    pathway = rownames(tt2),
    logFC = tt2$logFC,
    AveExpr = tt2$AveExpr,
    t = tt2$t,
    P.Value = tt2$P.Value,
    adj.P.Val = tt2$adj.P.Val,
    B = tt2$B,
    stringsAsFactors = FALSE
  )
  list(result=res2, mode="dupCorr")
}

# ---- Multi-contrast model (single clean fit per pathway DB) ----
run_multicontrast <- function(pbk, st, require_full_paired=TRUE, allow_fallback=TRUE) {
  
  stopifnot(all(c("subject","immune_status","treatment") %in% colnames(st)))
  
  st <- st %>%
    mutate(
      immune_status = factor(immune_status, levels=c("NonImmune","Immune")),
      treatment     = factor(treatment, levels=c("Untreated","Treated")),
      subject       = factor(subject)
    )
  
  if (!all(c("Untreated","Treated") %in% unique(as.character(st$treatment)))) {
    return(list(error="missing_treatment_levels(need_Untreated_and_Treated)"))
  }
  
  subj_full <- st %>%
    distinct(subject, immune_status, treatment) %>%
    count(subject) %>%
    filter(n >= 4) %>%
    pull(subject) %>%
    as.character()
  
  use_paired <- require_full_paired && length(subj_full) >= 2
  
  fit_with_contrasts <- function(pbk_use, st_use, use_subject_term, mode_label) {
    
    st_use <- st_use %>%
      mutate(
        subject = droplevels(subject),
        immune_status = droplevels(immune_status),
        treatment = droplevels(treatment)
      )
    
    if (use_subject_term) {
      design <- model.matrix(~ subject + immune_status * treatment, data=st_use)
    } else {
      design <- model.matrix(~ immune_status * treatment, data=st_use)
    }
    
    colnames(design) <- make.names(colnames(design))
    int_term <- "immune_statusImmune.treatmentTreated"
    
    needed <- c("immune_statusImmune", "treatmentTreated", int_term)
    missing_cols <- setdiff(needed, colnames(design))
    if (length(missing_cols) > 0) {
      return(list(error=paste0(
        "design_missing_terms: ", paste(missing_cols, collapse=", "),
        " | design_cols=", paste(colnames(design), collapse="|")
      )))
    }
    
    fit <- lmFit(pbk_use, design)
    fit <- tryCatch(eBayes(fit, trend=TRUE, robust=TRUE), error=function(e) NULL)
    if (is.null(fit)) return(list(error=paste0("eBayes_failed_", mode_label)))
    
    cm <- makeContrasts(
      Immune_vs_NonImmune_Untreated = immune_statusImmune,
      Immune_vs_NonImmune_Treated   = immune_statusImmune + `immune_statusImmune.treatmentTreated`,
      Immune_Treatment_Effect       = treatmentTreated + `immune_statusImmune.treatmentTreated`,
      NonImmune_Treatment_Effect    = treatmentTreated,
      Interaction_DeltaDelta        = `immune_statusImmune.treatmentTreated`,
      levels = design
    )
    
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- tryCatch(eBayes(fit2, trend=TRUE, robust=TRUE), error=function(e) NULL)
    if (is.null(fit2)) return(list(error=paste0("eBayes_failed_contrasts_", mode_label)))
    
    out <- lapply(colnames(cm), function(ct) {
      tt <- topTable(fit2, coef=ct, number=Inf, sort.by="none")
      data.frame(
        pathway = rownames(tt),
        logFC = tt$logFC,
        AveExpr = tt$AveExpr,
        t = tt$t,
        P.Value = tt$P.Value,
        adj.P.Val = tt$adj.P.Val,
        B = tt$B,
        contrast = ct,
        stringsAsFactors = FALSE
      )
    })
    
    list(results=bind_rows(out), mode=mode_label)
  }
  
  # Paired/subject model (preferred)
  if (use_paired) {
    st_use <- st %>% filter(as.character(subject) %in% subj_full)
    pbk_use <- pbk[, rownames(st_use), drop=FALSE]
    
    res <- fit_with_contrasts(pbk_use, st_use, use_subject_term=TRUE, mode_label="paired_subject")
    if (is.null(res$error)) return(res)
    
    if (!allow_fallback) return(res)
    logmsg("WARNING: paired subject model failed (%s). Falling back to dupCorr.", res$error)
  } else {
    if (!allow_fallback) {
      return(list(error=sprintf("insufficient_fully_paired_subjects(%d)", length(subj_full))))
    }
  }
  
  # Fallback: duplicateCorrelation (unbalanced subjects)
  st_use <- st %>%
    mutate(
      subject = droplevels(subject),
      immune_status = droplevels(immune_status),
      treatment = droplevels(treatment)
    )
  
  design <- model.matrix(~ immune_status * treatment, data=st_use)
  colnames(design) <- make.names(colnames(design))
  
  subject <- st_use$subject
  corfit <- duplicateCorrelation(pbk, design, block=subject)
  
  fit <- lmFit(pbk, design, block=subject, correlation=corfit$consensus)
  fit <- tryCatch(eBayes(fit, trend=TRUE, robust=TRUE), error=function(e) NULL)
  if (is.null(fit)) return(list(error="eBayes_failed_dupCorr"))
  
  int_term <- "immune_statusImmune.treatmentTreated"
  needed <- c("immune_statusImmune", "treatmentTreated", int_term)
  missing_cols <- setdiff(needed, colnames(design))
  if (length(missing_cols) > 0) {
    return(list(error=paste0(
      "design_missing_terms_dupCorr: ", paste(missing_cols, collapse=", "),
      " | design_cols=", paste(colnames(design), collapse="|")
    )))
  }
  
  cm <- makeContrasts(
    Immune_vs_NonImmune_Untreated = immune_statusImmune,
    Immune_vs_NonImmune_Treated   = immune_statusImmune + `immune_statusImmune.treatmentTreated`,
    Immune_Treatment_Effect       = treatmentTreated + `immune_statusImmune.treatmentTreated`,
    NonImmune_Treatment_Effect    = treatmentTreated,
    Interaction_DeltaDelta        = `immune_statusImmune.treatmentTreated`,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- tryCatch(eBayes(fit2, trend=TRUE, robust=TRUE), error=function(e) NULL)
  if (is.null(fit2)) return(list(error="eBayes_failed_contrasts_dupCorr"))
  
  out <- lapply(colnames(cm), function(ct) {
    tt <- topTable(fit2, coef=ct, number=Inf, sort.by="none")
    data.frame(
      pathway = rownames(tt),
      logFC = tt$logFC,
      AveExpr = tt$AveExpr,
      t = tt$t,
      P.Value = tt$P.Value,
      adj.P.Val = tt$adj.P.Val,
      B = tt$B,
      contrast = ct,
      stringsAsFactors = FALSE
    )
  })
  
  list(results=bind_rows(out), mode="dupCorr_unbalanced")
}

# =============================================================================
# STAGE 4: Analysis 1 — CD45+ vs CD45- per treatment
# =============================================================================
stage(4, "Analysis 1: CD45+ vs CD45- per treatment")

analyze_cd45_per_treatment <- function(escape_mat, pathway_type, meta_all) {
  treatment_levels <- sort(unique(meta_all[[TIME_COL]]))
  logmsg("Analyzing %s CD45 comparisons for treatments: %s", pathway_type, paste(treatment_levels, collapse=", "))
  
  all_results <- list()
  skip_rows <- list()
  
  for (trt in treatment_levels) {
    trt_safe <- sanitize_filename(trt)
    logmsg("  %s - %s", pathway_type, trt)
    
    meta_trt <- meta_all %>% filter(.data[[TIME_COL]] == trt)
    built <- build_pseudobulk(escape_mat, meta_trt)
    
    if (!is.null(built$error)) {
      skip_rows[[length(skip_rows)+1]] <- data.frame(
        pathway_type=pathway_type, analysis="CD45pos_vs_CD45neg", treatment=trt,
        status="SKIP", reason=built$error, stringsAsFactors=FALSE
      )
      next
    }
    
    pbk <- built$pbk
    st  <- built$st
    
    # QC
    write_build_qc(pathway_type, glue("CD45_{trt_safe}"), st, pbk)
    
    lim <- run_limma_cd45(pbk, st, allow_fallback = ALLOW_UNPAIRED_FALLBACK)
    if (!is.null(lim$error)) {
      skip_rows[[length(skip_rows)+1]] <- data.frame(
        pathway_type=pathway_type, analysis="CD45pos_vs_CD45neg", treatment=trt,
        status="SKIP", reason=lim$error, stringsAsFactors=FALSE
      )
      next
    }
    
    df <- lim$result %>% mutate(
      pathway_type=pathway_type,
      treatment=trt,
      method=paste0("limma_", lim$mode),
      analysis="CD45pos_vs_CD45neg"
    )
    
    all_results[[paste0(pathway_type, "_", trt)]] <- df
    
    out_csv <- file.path(DIR_DEG, glue("{pathway_type}_CD45comparison_{trt_safe}.csv"))
    fwrite(df, out_csv)
    logmsg("Wrote: %s", out_csv)
    
    out_png <- file.path(DIR_PLOT, glue("{pathway_type}_volcano_CD45_{trt_safe}.png"))
    plot_volcano(df, glue("{pathway_type}: CD45+ vs CD45- ({trt})"), out_png)
    
    skip_rows[[length(skip_rows)+1]] <- data.frame(
      pathway_type=pathway_type, analysis="CD45pos_vs_CD45neg", treatment=trt,
      status="OK", reason="", stringsAsFactors=FALSE
    )
  }
  
  list(results=all_results, skips=skip_rows)
}

reactome_cd45 <- analyze_cd45_per_treatment(escape_reactome, "Reactome", meta_all)

gobp_cd45 <- list(results=list(), skips=list())
if (!is.null(escape_gobp)) {
  gobp_cd45 <- analyze_cd45_per_treatment(escape_gobp, "GOBP", meta_all)
}

# =============================================================================
# STAGE 5: Analysis 2 — Full multi-contrast model (single fit per pathway DB)
# =============================================================================
stage(5, "Analysis 2: Full multi-contrast model")

analyze_multicontrast <- function(escape_mat, pathway_type, meta_all) {
  logmsg("Analyzing %s multi-contrast model", pathway_type)
  
  built <- build_pseudobulk(escape_mat, meta_all)
  if (!is.null(built$error)) return(list(error=built$error))
  
  pbk <- built$pbk
  st  <- built$st
  
  # QC
  write_build_qc(pathway_type, "MultiContrast_ALL", st, pbk)
  
  mc <- run_multicontrast(
    pbk, st,
    require_full_paired = TRUE,
    allow_fallback      = ALLOW_UNPAIRED_FALLBACK
  )
  if (!is.null(mc$error)) return(list(error=mc$error))
  
  res_all <- mc$results %>%
    mutate(
      pathway_type = pathway_type,
      method = paste0("limma_", mc$mode),
      analysis = "MultiContrast_subject_plus_interaction"
    )
  
  out_csv <- file.path(DIR_DEG, glue("{pathway_type}_AllContrasts.csv"))
  fwrite(res_all, out_csv)
  logmsg("Wrote: %s", out_csv)
  
  for (ct in unique(res_all$contrast)) {
    df_ct <- res_all %>% filter(contrast == ct)
    out_png <- file.path(DIR_PLOT, glue("{pathway_type}_volcano_{sanitize_filename(ct)}.png"))
    plot_volcano(df_ct, glue("{pathway_type}: {ct}"), out_png)
  }
  
  list(result=res_all)
}

reactome_mc <- analyze_multicontrast(escape_reactome, "Reactome", meta_all)

gobp_mc <- list(result=NULL)
if (!is.null(escape_gobp)) {
  gobp_mc <- analyze_multicontrast(escape_gobp, "GOBP", meta_all)
}

# =============================================================================
# STAGE 6: Save combined reports + summaries
# =============================================================================
stage(6, "Save reports")

# ---- Combined CD45 comparison results ----
all_cd45_results <- c(reactome_cd45$results, gobp_cd45$results)
if (length(all_cd45_results) > 0) {
  combined_cd45 <- bind_rows(all_cd45_results)
  out_cd45 <- file.path(DIR_DEG, "ALL_CD45pos_vs_CD45neg.csv")
  fwrite(combined_cd45, out_cd45)
  logmsg("Wrote: %s", out_cd45)
  
  top_cd45 <- combined_cd45 %>%
    filter(adj.P.Val <= FDR_CUTOFF, abs(logFC) >= LOGFC_CUTOFF) %>%
    arrange(pathway_type, treatment, adj.P.Val)
  
  out_top_cd45 <- file.path(DIR_DEG, "TOP_CD45_significant.csv")
  fwrite(top_cd45, out_top_cd45)
  logmsg("Wrote: %s", out_top_cd45)
  
  # Summary counts
  sum_cd45 <- combined_cd45 %>%
    mutate(sig = (adj.P.Val <= FDR_CUTOFF & abs(logFC) >= LOGFC_CUTOFF)) %>%
    group_by(pathway_type, treatment, method, analysis) %>%
    summarise(
      n_pathways = dplyr::n(),
      n_sig = sum(sig, na.rm=TRUE),
      .groups="drop"
    ) %>%
    arrange(pathway_type, treatment, method)
  out_sum_cd45 <- file.path(DIR_DEG, "SUMMARY_CD45_counts.csv")
  fwrite(sum_cd45, out_sum_cd45)
  logmsg("Wrote: %s", out_sum_cd45)
}

# ---- Combined multi-contrast results ----
mc_list <- list()
if (!is.null(reactome_mc$result)) mc_list <- c(mc_list, list(reactome_mc$result))
if (!is.null(gobp_mc$result))     mc_list <- c(mc_list, list(gobp_mc$result))

if (length(mc_list) > 0) {
  combined_mc <- bind_rows(mc_list)
  out_mc <- file.path(DIR_DEG, "ALL_AllContrasts.csv")
  fwrite(combined_mc, out_mc)
  logmsg("Wrote: %s", out_mc)
  
  for (ct in unique(combined_mc$contrast)) {
    top_ct <- combined_mc %>%
      filter(contrast == ct, adj.P.Val <= FDR_CUTOFF, abs(logFC) >= LOGFC_CUTOFF) %>%
      arrange(pathway_type, adj.P.Val)
    
    out_top_ct <- file.path(DIR_DEG, glue("TOP_{sanitize_filename(ct)}.csv"))
    fwrite(top_ct, out_top_ct)
    logmsg("Wrote: %s", out_top_ct)
  }
  
  # Summary counts per contrast
  sum_mc <- combined_mc %>%
    mutate(sig = (adj.P.Val <= FDR_CUTOFF & abs(logFC) >= LOGFC_CUTOFF)) %>%
    group_by(pathway_type, contrast, method, analysis) %>%
    summarise(
      n_pathways = dplyr::n(),
      n_sig = sum(sig, na.rm=TRUE),
      .groups="drop"
    ) %>%
    arrange(pathway_type, contrast, method)
  
  out_sum_mc <- file.path(DIR_DEG, "SUMMARY_MultiContrast_counts.csv")
  fwrite(sum_mc, out_sum_mc)
  logmsg("Wrote: %s", out_sum_mc)
  
} else {
  logmsg("No multi-contrast results to combine.")
  if (!is.null(reactome_mc$error)) logmsg("Reactome multi-contrast error: %s", reactome_mc$error)
  if (!is.null(gobp_mc$error))     logmsg("GOBP multi-contrast error: %s", gobp_mc$error)
}

# ---- Skip report ----
all_skips <- bind_rows(c(reactome_cd45$skips, gobp_cd45$skips))
if (nrow(all_skips) > 0) {
  out_skip <- file.path(OUT_ROOT, "SKIP_REPORT.csv")
  fwrite(all_skips, out_skip)
  logmsg("Wrote: %s", out_skip)
}

stage(7, "Done")
logmsg("Output: %s", OUT_ROOT)
logmsg("Done!")
