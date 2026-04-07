#!/usr/bin/env Rscript
################################################################################
# Comprehensive Single-Cell RNA-seq Analysis Pipeline – PBMC VERSION
# Experiment 301-4  (begun Dec 17th)
#
# Design:
#   16 PBMC samples from 12 subjects loaded across 6 lanes (all lanes identical)
#   Pool 1 (HTO_3): 4 APECED patients T1 + 4 HC + 4 Severe
#   Pool 2 (HTO_4): 4 APECED patients T2 + 4 HC + 4 Severe
#   HTO hashtags (Antibody Capture inside .h5) separate the two pools.
#   HTODemux() is run PER-LANE before merging; pool identity comes from hash.ID.
#   TimePoint is derived from HTO pool membership – NOT from lane number.
#
# Outputs: Raw object + Final annotated object (2 RDS files)
#
# PBMC CHANGE markers: every modification from the ileum script is flagged
# with a comment beginning  ## PBMC CHANGE:
################################################################################

## PBMC CHANGE: hashtag-to-pool reference (kept at top for easy editing)
#   HTO_3  →  Pool 1  →  T1 for APECED patients
#   HTO_4  →  Pool 2  →  T2 for APECED patients
#   HC and Severe are present in both pools with no timepoint distinction

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(cowplot)
  library(readr)
  library(stringr)
  library(purrr)
  library(patchwork)
  library(scales)
  library(Matrix)
  library(data.table)
  library(SingleR)
  library(celldex, lib = "../../lib")
  library(BiocParallel)
  library(RColorBrewer)
})

set.seed(1234)

################################################################################
# CONFIGURATION
################################################################################

## PBMC CHANGE: project label updated
PROJECT  <- "APECED_PBMC"

BASE_DIR <- "/gpfs/gsfs12/users/cheungf/ACEPD/PBMC"
OUT_DIR  <- file.path(BASE_DIR, "PBMC_PIPELINE_OUTPUT")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## PBMC CHANGE: subject vectors now include HC and Severe cohorts
# APECED patients who have T1 (Pool1/HTO_3) and T2 (Pool2/HTO_4)
DONORS <- c("P0002345", "P0001171", "P0001358", "P0001352")

# Severe patients – present in both pools, labelled Severe_Pool1 / Severe_Pool2
SEVERE <- c("P0001140", "P0011505", "P0005819", "P0014164")

# Healthy controls – spike-in in both pools, collapsed to a single label
# CTRL_RAW: the individual IDs as they appear in SNP assignment files
CTRL_RAW_IDS <- c("HC.01", "HC.02", "HC.03", "HC.04")
# CTRL_LAB: harmonised label used throughout the script
CTRL_LAB     <- "HC"

## PBMC CHANGE: HTO pool identifiers as written by HTODemux into hash.ID
HTO_POOL1 <- "HTO_3"   # Pool 1 → T1 for APECED patients
HTO_POOL2 <- "HTO_4"   # Pool 2 → T2 for APECED patients

QC_THRESHOLDS <- list(
  nCount_min       = 1200,
  nCount_max       = 50000,
  nFeature_min     = 300,
  nFeature_max     = 7000,
  percent_mito_max = 60,
  percent_ribo_max = 25,
  percent_hb_max   = 0.1
)

ANNOTATION_CONFIG <- list(
  p_adj_max        = 0.05,
  only_pos         = TRUE,
  min_pct          = 0.25,
  logfc_threshold  = 0.25,
  tolerance        = 0.10,
  force_recompute  = FALSE
)

################################################################################
# UTILITY FUNCTIONS  (unchanged from ileum script)
################################################################################

log_status <- function(stage, message, level = "INFO") {
  icons <- list("INFO" = "ℹ️", "SUCCESS" = "✅", "WARNING" = "⚠️",
                "ERROR" = "❌", "PROGRESS" = "▶️")
  icon      <- icons[[level]]
  timestamp <- format(Sys.time(), "%H:%M:%S")
  cat(sprintf("[%s] %s [%s] %s\n", timestamp, icon, stage, message))
}

section_header <- function(title) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("  %s\n", title))
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
}

save_plot_quietly <- function(filename, plot, ...) {
  suppressMessages(ggsave(filename, plot, ...))
  log_status("PLOT", sprintf("Saved: %s", basename(filename)), "SUCCESS")
}

check_and_join_layers <- function(obj, stage_name) {
  current_layers <- SeuratObject::Layers(obj[["RNA"]])
  data_layers    <- grep("^data",   current_layers, value = TRUE)
  counts_layers  <- grep("^counts", current_layers, value = TRUE)
  needs_join     <- (length(data_layers) > 1) || (length(counts_layers) > 1)
  if (needs_join) {
    log_status(stage_name, "Multiple layers detected - joining...", "INFO")
    if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
      obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
      log_status(stage_name, "Layers joined successfully", "SUCCESS")
    }
  }
  return(obj)
}

join_if_multi <- function(obj, assay = "RNA", stage_name = "JOIN") {
  if (!("JoinLayers" %in% getNamespaceExports("SeuratObject"))) return(obj)
  layers        <- tryCatch(SeuratObject::Layers(obj[[assay]]),
                            error = function(e) character(0))
  data_layers   <- grep("^data",   layers, value = TRUE)
  counts_layers <- grep("^counts", layers, value = TRUE)
  if (length(data_layers) > 1 || length(counts_layers) > 1) {
    log_status(stage_name, "Collapsing layers...", "PROGRESS")
    obj <- SeuratObject::JoinLayers(obj, assay = assay)
    log_status(stage_name, "Layers joined", "SUCCESS")
  }
  obj
}

################################################################################
# STAGE 1: DATA LOADING & INTEGRATION
## PBMC CHANGE: reads both Gene Expression AND Antibody Capture from each h5,
##   creates an HTO assay, normalises it with CLR, and runs HTODemux() per-lane
##   before merging.  The GEX assay handling is identical to the ileum script.
################################################################################

stage1_load_and_integrate <- function() {
  section_header("STAGE 1: DATA LOADING & INTEGRATION (with per-lane HTODemux)")
  
  h5_paths <- c(
    "23-301-4_GEX1" = "cellranger_outs/filtered_feature_bc_matrix_23-301-4_GEX1.h5",
    "23-301-4_GEX2" = "cellranger_outs/filtered_feature_bc_matrix_23-301-4_GEX2.h5",
    "23-301-4_GEX3" = "cellranger_outs/filtered_feature_bc_matrix_23-301-4_GEX3.h5",
    "23-301-4_GEX4" = "cellranger_outs/filtered_feature_bc_matrix_23-301-4_GEX4.h5",
    "23-301-4_GEX5" = "cellranger_outs/filtered_feature_bc_matrix_23-301-4_GEX5.h5",
    "23-301-4_GEX6" = "cellranger_outs/filtered_feature_bc_matrix_23-301-4_GEX6.h5"
  )
  
  log_status("LOAD", sprintf("Loading %d GEX samples with HTO...", length(h5_paths)),
             "PROGRESS")
  
  objs <- list()
  
  for (s in names(h5_paths)) {
    p <- h5_paths[[s]]
    if (!file.exists(p)) {
      log_status("LOAD", sprintf("Missing file: %s", p), "ERROR")
      stop("Missing required file")
    }
    
    log_status("LOAD", sprintf("Reading %s", s), "INFO")
    m <- Read10X_h5(p, use.names = TRUE)
    
    # ── Gene Expression ──────────────────────────────────────────────────────
    ge <- if (is.list(m) && "Gene Expression" %in% names(m)) {
      m[["Gene Expression"]]
    } else if (inherits(m, "dgCMatrix")) {
      m
    } else {
      stop("Unrecognised Read10X_h5 structure for GEX in ", s)
    }
    
    clean_bc    <- sub("-.*$", "", colnames(ge))
    colnames(ge) <- paste0(clean_bc, "-", s)
    
    sobj <- CreateSeuratObject(counts = ge, project = PROJECT)
    sobj$orig.ident <- s
    sobj <- join_if_multi(sobj, stage_name = paste0("JOIN_GEX_", s))
    
    ## PBMC CHANGE: extract Antibody Capture matrix and add as HTO assay -----
    if (is.list(m) && "Antibody Capture" %in% names(m)) {
      hto_mat <- m[["Antibody Capture"]]
      
      # align barcodes to the same suffix used for GEX
      hto_clean_bc    <- sub("-.*$", "", colnames(hto_mat))
      colnames(hto_mat) <- paste0(hto_clean_bc, "-", s)
      
      # keep only barcodes present in the GEX object
      shared_bc <- intersect(colnames(sobj), colnames(hto_mat))
      if (length(shared_bc) < ncol(sobj)) {
        log_status("HTO",
                   sprintf("%s: %d/%d barcodes have HTO data",
                           s, length(shared_bc), ncol(sobj)), "WARNING")
      }
      
      sobj <- subset(sobj, cells = shared_bc)
      hto_mat <- hto_mat[, shared_bc, drop = FALSE]
      
      sobj[["HTO"]] <- CreateAssayObject(counts = hto_mat)
      
      # CLR normalisation across cells (margin = 2) as recommended by Seurat
      sobj <- NormalizeData(sobj, assay = "HTO",
                            normalization.method = "CLR", margin = 2,
                            verbose = FALSE)
      
      # HTODemux per-lane: positive.quantile can be tuned if needed
      log_status("HTO", sprintf("Running HTODemux on %s ...", s), "PROGRESS")
      sobj <- tryCatch(
        HTODemux(sobj, assay = "HTO", positive.quantile = 0.99, verbose = FALSE),
        error = function(e) {
          log_status("HTO",
                     sprintf("HTODemux failed for %s: %s", s, e$message),
                     "WARNING")
          # add NA placeholder columns so downstream code never breaks
          sobj$HTO_maxID                  <- NA_character_
          sobj$HTO_secondID               <- NA_character_
          sobj$HTO_margin                 <- NA_real_
          sobj$HTO_classification         <- NA_character_
          sobj$HTO_classification.global  <- NA_character_
          sobj$hash.ID                    <- NA_character_
          sobj
        }
      )
      
      n_pool1 <- sum(sobj$hash.ID == HTO_POOL1, na.rm = TRUE)
      n_pool2 <- sum(sobj$hash.ID == HTO_POOL2, na.rm = TRUE)
      n_doub  <- sum(sobj$HTO_classification.global == "Doublet", na.rm = TRUE)
      n_neg   <- sum(sobj$HTO_classification.global == "Negative", na.rm = TRUE)
      log_status("HTO",
                 sprintf("%s → Pool1(%s):%d  Pool2(%s):%d  Doublet:%d  Negative:%d",
                         s, HTO_POOL1, n_pool1, HTO_POOL2, n_pool2, n_doub, n_neg),
                 "SUCCESS")
    } else {
      log_status("HTO",
                 sprintf("%s: no Antibody Capture feature found – HTO skipped", s),
                 "WARNING")
      sobj$hash.ID                   <- NA_character_
      sobj$HTO_classification.global <- NA_character_
    }
    ## END PBMC CHANGE --------------------------------------------------------
    
    log_status("LOAD", sprintf("%s: %d cells", s, ncol(sobj)), "SUCCESS")
    objs[[s]] <- sobj
  }
  
  log_status("MERGE", "Merging all samples...", "PROGRESS")
  allsamples <- Reduce(function(a, b) merge(a, b, project = PROJECT), objs)
  allsamples <- join_if_multi(allsamples, stage_name = "JOIN_MERGED")
  
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    allsamples <- SeuratObject::JoinLayers(allsamples)
    log_status("MERGE", "Joined layers (Seurat v5)", "SUCCESS")
  }
  
  log_status("MERGE",
             sprintf("Total cells after per-lane HTODemux merge: %s",
                     format(ncol(allsamples), big.mark = ",")),
             "SUCCESS")
  
  return(allsamples)
}

################################################################################
# STAGE 2: METADATA INTEGRATION (WITH MULTI-TOOL DOUBLETS)
## PBMC CHANGE: HTO pool metadata section added (uses hash.ID already in the
##   object from Stage 1).  SNP demux, ambient RNA, and multi-tool doublet
##   sections are unchanged from the ileum script except for the file glob
##   pattern which now targets 23-301-4_* files.
################################################################################

stage2_add_metadata <- function(obj) {
  section_header("STAGE 2: METADATA INTEGRATION")
  
  # ---------------------------
  # QC metrics
  # ---------------------------
  log_status("QC", "Calculating QC metrics...", "PROGRESS")
  obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$percent.ribo <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  hb_genes         <- grep("^HB", rownames(obj), value = TRUE)
  hb_genes         <- hb_genes[!grepl("^HBP", hb_genes)]
  obj$percent.hb   <- PercentageFeatureSet(obj, features = hb_genes)
  log_status("QC", "QC metrics calculated", "SUCCESS")
  
  # ---------------------------
  ## PBMC CHANGE: HTO pool metadata
  # hash.ID is already present from HTODemux() in Stage 1.
  # Derive a clean Pool column and flag HTO doublets / negatives.
  # ---------------------------
  log_status("HTO_META", "Deriving pool metadata from HTODemux results...", "PROGRESS")
  
  if ("hash.ID" %in% colnames(obj@meta.data)) {
    obj$HTO_pool <- dplyr::case_when(
      obj$hash.ID == HTO_POOL1                         ~ "Pool1",
      obj$hash.ID == HTO_POOL2                         ~ "Pool2",
      obj$HTO_classification.global == "Doublet"       ~ "HTO_Doublet",
      obj$HTO_classification.global == "Negative"      ~ "HTO_Negative",
      TRUE                                              ~ "HTO_Unknown"
    )
    
    pool_tbl <- table(obj$HTO_pool)
    for (nm in names(pool_tbl)) {
      log_status("HTO_META",
                 sprintf("  %s: %s cells", nm,
                         format(pool_tbl[[nm]], big.mark = ",")), "INFO")
    }
    log_status("HTO_META", "HTO_pool column added", "SUCCESS")
  } else {
    log_status("HTO_META",
               "hash.ID not found – HTO pool assignment skipped", "WARNING")
    obj$HTO_pool <- NA_character_
  }
  ## END PBMC CHANGE ----------------------------------------------------------
  
  # ---------------------------
  # SNP demux
  # ---------------------------
  log_status("SNP", "Loading SNP demultiplexing data...", "PROGRESS")
  ## PBMC CHANGE: glob updated to 23-301-4 experiment ID
  assign_files <- Sys.glob("./DEMUX/*/assignments_refined.tsv")
  
  if (length(assign_files) > 0) {
    read_one <- function(p) {
      ## PBMC CHANGE: regex updated to match 23-301-4_GEX* lane names
      smp     <- stringr::str_match(p, "(23-301-4_GEX\\d+)_demux")[, 2]
      df      <- readr::read_tsv(p,
                                 col_types  = readr::cols(.default = readr::col_character()),
                                 show_col_types = FALSE)
      subj_col <- setdiff(names(df), "BARCODE")[1]
      df %>%
        dplyr::transmute(
          core_bc = sub("-.*$", "", BARCODE),
          subject = .data[[subj_col]],
          cell    = paste0(core_bc, "-", smp)
        ) %>%
        dplyr::select(cell, subject)
    }
    
    assign_tbl <- purrr::map_dfr(assign_files, read_one) %>%
      dplyr::distinct(cell, .keep_all = TRUE)
    
    meta <- obj@meta.data %>% dplyr::mutate(cell = rownames(.))
    
    meta2 <- meta %>%
      dplyr::left_join(assign_tbl, by = "cell", suffix = c("", "_snp")) %>%
      dplyr::mutate(
        Lane     = stringr::str_match(orig.ident, "GEX(\\d+)")[, 2],
        SNP_call = dplyr::case_when(
          is.na(subject)                                  ~ "Unassigned",
          stringr::str_detect(as.character(subject), "\\+") ~ "Doublet",
          TRUE                                            ~ "Singleton"
        )
      )
    meta2 <- meta2[match(colnames(obj), meta2$cell), , drop = FALSE]
    
    obj$subject  <- meta2$subject
    obj$Lane     <- meta2$Lane
    obj$SNP_call <- meta2$SNP_call
    log_status("SNP",
               sprintf("Assigned: %d cells", sum(!is.na(obj$subject))),
               "SUCCESS")
  } else {
    log_status("SNP", "No SNP assignment files found; skipping.", "WARN")
  }
  
  # Derive Lane from orig.ident regardless of whether SNP demux ran.
  # orig.ident is always "23-301-4_GEXn" so we just pull the trailing digit(s).
  if (!"Lane" %in% colnames(obj@meta.data)) {
    obj$Lane <- stringr::str_match(as.character(obj$orig.ident), "GEX(\\d+)")[, 2]
    log_status("SNP",
               sprintf("Lane derived from orig.ident: %s",
                       paste(sort(unique(obj$Lane)), collapse = ", ")), "INFO")
  }
  
  # ---------------------------
  # Ambient RNA
  # ---------------------------
  log_status("AMBIENT", "Loading ambient RNA scores...", "PROGRESS")
  ## PBMC CHANGE: glob updated to 23-301-4 experiment ID
  ambient_files <- Sys.glob("AMBIENT/*GEX*_ambient_scores.tsv")
  
  if (length(ambient_files) > 0) {
    read_ambient <- function(f) {
      ## PBMC CHANGE: regex updated to match 23-301-4_GEX* lane names
      smp <- stringr::str_match(basename(f), "(23-301-4_GEX\\d+)_ambient")[, 2]
      df  <- readr::read_tsv(f,
                             col_types = readr::cols(.default = readr::col_character()),
                             show_col_types = FALSE)
      df %>%
        dplyr::mutate(
          core_bc             = sub("-.*$", "", barcode),
          cell                = paste0(core_bc, "-", smp),
          ambient_rna_score   = suppressWarnings(as.numeric(ambient_rna_score)),
          percent_top_ambient = suppressWarnings(as.numeric(percent_top_ambient)),
          library_complexity  = suppressWarnings(as.numeric(library_complexity))
        ) %>%
        dplyr::select(cell, ambient_rna_score, percent_top_ambient, library_complexity)
    }
    
    ambient_tbl <- purrr::map_dfr(ambient_files, read_ambient) %>%
      dplyr::distinct(cell, .keep_all = TRUE)
    
    meta3 <- obj@meta.data %>%
      dplyr::mutate(cell = rownames(.)) %>%
      dplyr::left_join(ambient_tbl, by = "cell", suffix = c("", ".ambient")) %>%
      { .[match(colnames(obj), .$cell), , drop = FALSE] }
    
    get_col <- function(base) {
      a <- paste0(base, ".ambient")
      if (a %in% names(meta3)) {
        if (base %in% names(meta3)) dplyr::coalesce(meta3[[a]], meta3[[base]])
        else meta3[[a]]
      } else if (base %in% names(meta3)) {
        meta3[[base]]
      } else {
        rep(NA_real_, ncol(obj))
      }
    }
    
    obj$ambient_rna_score   <- get_col("ambient_rna_score")
    obj$percent_top_ambient <- get_col("percent_top_ambient")
    obj$library_complexity  <- get_col("library_complexity")
    
    log_status("AMBIENT",
               sprintf("Added scores for %d cells",
                       sum(!is.na(obj$ambient_rna_score))),
               "SUCCESS")
  } else {
    log_status("AMBIENT", "No ambient score files found; skipping.", "WARN")
  }
  
  # ---------------------------
  # MULTI-TOOL DOUBLET DETECTION
  # ---------------------------
  log_status("DOUBLET", "Loading multi-tool doublet detection results...", "PROGRESS")
  
  ## PBMC CHANGE: glob updated to 23-301-4 experiment ID
  doublet_files <- Sys.glob("./DOUBLET/23-301-4_*_combined_DemuxOnly.tsv")
  if (length(doublet_files) == 0)
    doublet_files <- Sys.glob("./DOUBLET/*combined_DemuxOnly.tsv")
  if (length(doublet_files) == 0) {
    fallback_file <- "/mnt/data/23-301-4_ALL_combined_DemuxOnly.tsv"
    if (file.exists(fallback_file)) doublet_files <- fallback_file
  }
  
  if (length(doublet_files) > 0) {
    
    read_doublets_multi <- function(f) {
      log_status("DOUBLET", sprintf("  Reading %s", basename(f)), "INFO")
      df <- readr::read_tsv(f,
                            col_types      = readr::cols(.default = readr::col_character()),
                            show_col_types = FALSE,
                            guess_max      = 1e6,
                            trim_ws        = TRUE)
      
      if (!"Barcode" %in% names(df))
        stop("Doublet TSV missing 'Barcode' column: ", f)
      
      df %>%
        dplyr::mutate(
          Barcode  = trimws(Barcode),
          core_bc  = toupper(sub("-.*$", "", Barcode)),
          Demuxalot_Individual_Assignment =
            if ("Demuxalot_Individual_Assignment" %in% names(.))
              trimws(Demuxalot_Individual_Assignment) else NA_character_,
          Demuxalot_DropletType =
            if ("Demuxalot_DropletType" %in% names(.))
              tolower(trimws(Demuxalot_DropletType)) else NA_character_,
          DoubletDetection_DropletType =
            if ("DoubletDetection_DropletType" %in% names(.))
              tolower(trimws(DoubletDetection_DropletType)) else NA_character_,
          scDblFinder_DropletType =
            if ("scDblFinder_DropletType" %in% names(.))
              tolower(trimws(scDblFinder_DropletType)) else NA_character_,
          scDblFinder_Score =
            if ("scDblFinder_Score" %in% names(.))
              suppressWarnings(as.numeric(scDblFinder_Score)) else NA_real_,
          scds_score =
            if ("scds_score" %in% names(.))
              suppressWarnings(as.numeric(scds_score)) else NA_real_,
          scds_DropletType =
            if ("scds_DropletType" %in% names(.))
              tolower(trimws(scds_DropletType)) else NA_character_,
          scrublet_DropletType =
            if ("scrublet_DropletType" %in% names(.))
              tolower(trimws(scrublet_DropletType)) else NA_character_,
          scrublet_Scores =
            if ("scrublet_Scores" %in% names(.))
              suppressWarnings(as.numeric(scrublet_Scores)) else NA_real_
        ) %>%
        dplyr::transmute(
          core_bc,
          Demuxalot_Individual     = Demuxalot_Individual_Assignment,
          Demuxalot_Droplet        = Demuxalot_DropletType,
          DoubletDetection_Droplet = DoubletDetection_DropletType,
          scDblFinder_Droplet      = scDblFinder_DropletType,
          scDblFinder_Score,
          scds_Score               = scds_score,
          scds_Droplet             = scds_DropletType,
          scrublet_Droplet         = scrublet_DropletType,
          scrublet_Score           = scrublet_Scores
        )
    }
    
    dbl_raw <- purrr::map_dfr(doublet_files, read_doublets_multi) %>%
      dplyr::mutate(core_bc = trimws(core_bc)) %>%
      dplyr::distinct(core_bc, .keep_all = TRUE)
    
    obj_names_uc  <- toupper(trimws(colnames(obj)))
    lanes_in_obj  <- unique(sub("^[^-]+-", "", obj_names_uc))
    
    cand_mat <- if (length(lanes_in_obj) > 0 && nrow(dbl_raw) > 0) {
      sapply(lanes_in_obj,
             function(ln) paste0(dbl_raw$core_bc, "-", ln) %in% obj_names_uc)
    } else {
      matrix(FALSE, nrow = nrow(dbl_raw), ncol = 0)
    }
    
    any_hit     <- if (ncol(cand_mat) > 0) rowSums(cand_mat) > 0 else
      rep(FALSE, nrow(dbl_raw))
    chosen_idx  <- if (ncol(cand_mat) > 0)
      max.col(cand_mat, ties.method = "first") else
        integer(nrow(dbl_raw))
    chosen_lane <- if (length(lanes_in_obj) > 0) lanes_in_obj[chosen_idx] else
      rep(NA_character_, nrow(dbl_raw))
    chosen_lane[!any_hit] <- NA_character_
    
    placed    <- any_hit
    n_matched <- sum(placed)
    
    if (n_matched == 0)
      log_status("DOUBLET",
                 "Matched 0 cells – barcodes may not overlap this object.", "WARN")
    
    dbl_tbl <- dbl_raw[placed, , drop = FALSE] %>%
      dplyr::mutate(cell = paste0(core_bc, "-", chosen_lane[placed])) %>%
      dplyr::select(-core_bc)
    
    droplet_cols <- intersect(
      c("Demuxalot_Droplet", "DoubletDetection_Droplet",
        "scDblFinder_Droplet", "scds_Droplet", "scrublet_Droplet"),
      colnames(dbl_tbl)
    )
    
    if (length(droplet_cols) > 0) {
      dbl_mat    <- as.matrix(dbl_tbl[droplet_cols] == "doublet")
      called_mat <- as.matrix(!is.na(dbl_tbl[droplet_cols]))
      n_doublet_calls <- rowSums(dbl_mat,    na.rm = TRUE)
      n_tools_called  <- rowSums(called_mat, na.rm = TRUE)
    } else {
      n_doublet_calls <- integer(nrow(dbl_tbl))
      n_tools_called  <- integer(nrow(dbl_tbl))
    }
    
    dbl_tbl <- dbl_tbl %>%
      dplyr::mutate(
        n_doublet_calls   = n_doublet_calls,
        n_tools_called    = n_tools_called,
        doublet_consensus = dplyr::case_when(
          n_tools_called == 0 ~ "unknown",
          n_doublet_calls / pmax(n_tools_called, 1L) >= 0.5 ~ "doublet",
          TRUE ~ "singlet"
        ),
        doublet_confidence = dplyr::if_else(
          n_tools_called > 0,
          pmax(n_doublet_calls, n_tools_called - n_doublet_calls) / n_tools_called,
          as.numeric(NA)
        )
      )
    
    meta_doublet <- obj@meta.data %>%
      dplyr::mutate(cell = rownames(.)) %>%
      dplyr::left_join(dbl_tbl, by = "cell", suffix = c("", ".dbl")) %>%
      { .[match(colnames(obj), .$cell), , drop = FALSE] }
    
    get_dbl <- function(base) {
      new <- paste0(base, ".dbl")
      if (new %in% names(meta_doublet) || base %in% names(meta_doublet))
        dplyr::coalesce(meta_doublet[[new]], meta_doublet[[base]])
      else
        rep(NA, ncol(obj))
    }
    
    obj$Demuxalot_Individual     <- get_dbl("Demuxalot_Individual")
    obj$Demuxalot_Droplet        <- get_dbl("Demuxalot_Droplet")
    obj$DoubletDetection_Droplet <- get_dbl("DoubletDetection_Droplet")
    obj$scDblFinder_Droplet      <- get_dbl("scDblFinder_Droplet")
    obj$scDblFinder_Score        <- get_dbl("scDblFinder_Score")
    obj$scds_Score               <- get_dbl("scds_Score")
    obj$scds_Droplet             <- get_dbl("scds_Droplet")
    obj$scrublet_Droplet         <- get_dbl("scrublet_Droplet")
    obj$scrublet_Score           <- get_dbl("scrublet_Score")
    obj$doublet_consensus        <- get_dbl("doublet_consensus")
    obj$doublet_confidence       <- get_dbl("doublet_confidence")
    obj$n_doublet_calls          <- get_dbl("n_doublet_calls")
    obj$n_tools_called           <- get_dbl("n_tools_called")
    
    log_status(
      "DOUBLET",
      sprintf("Attached doublet calls to %s / %s cells (cores placed: %s)",
              format(sum(!is.na(obj$Demuxalot_Droplet) |
                           !is.na(obj$DoubletDetection_Droplet) |
                           !is.na(obj$scDblFinder_Droplet)      |
                           !is.na(obj$scds_Droplet)             |
                           !is.na(obj$scrublet_Droplet)),
                     big.mark = ","),
              format(ncol(obj), big.mark = ","),
              format(n_matched, big.mark = ",")),
      if (n_matched > 0) "SUCCESS" else "WARN"
    )
    
    cat("\n")
    for (tool_col in c("DoubletDetection_Droplet", "scDblFinder_Droplet",
                       "scds_Droplet", "scrublet_Droplet")) {
      if (tool_col %in% colnames(obj@meta.data)) {
        tool_name <- sub("_Droplet$", "", tool_col)
        denom     <- sum(!is.na(obj@meta.data[[tool_col]]))
        if (denom > 0) {
          n_dbl   <- sum(obj@meta.data[[tool_col]] == "doublet", na.rm = TRUE)
          pct_dbl <- 100 * n_dbl / denom
          log_status("DOUBLET",
                     sprintf("  %s: %.1f%%", tool_name, pct_dbl), "INFO")
        } else {
          log_status("DOUBLET", sprintf("  %s: no calls", tool_name), "INFO")
        }
      }
    }
    cat("\n")
    
  } else {
    log_status("DOUBLET", "No doublet files found; skipping.", "WARN")
  }
  
  return(obj)
}

################################################################################
# STAGE 3: SINGLETON FILTERING & HARMONIZATION
## PBMC CHANGE (entire function rewritten):
##   1. Whitelist now includes DONORS, SEVERE, and CTRL_LAB
##   2. Subject_fixed harmonises the four individual HC IDs → "HC"
##      and keeps Severe IDs as-is (TimePoint logic below separates them)
##   3. TimePoint_fixed is derived from HTO_pool (Pool1 / Pool2), NOT lane
##      - APECED donors: Pool1 → "Before", Pool2 → "After"
##      - HC:            always "Spike" (same as ileum)
##      - Severe:        Pool1 → "Severe_Pool1", Pool2 → "Severe_Pool2"
##   4. Treatment gains "Severe_Pool1", "Severe_Pool2", and "HC" arms
################################################################################

stage3_filter_and_harmonize <- function(obj) {
  section_header("STAGE 3: SINGLETON FILTERING & HARMONIZATION")
  
  md <- obj@meta.data
  md$Subject_raw   <- as.character(md$subject)
  md$Subject_fixed <- md$Subject_raw
  
  # ── Harmonise subject labels ──────────────────────────────────────────────
  # Individual HC IDs → single label "HC"
  md$Subject_fixed[md$Subject_fixed %in% CTRL_RAW_IDS] <- CTRL_LAB
  # Catch any alternative raw string the SNP caller may emit for controls
  md$Subject_fixed[md$Subject_fixed == "Control"]       <- CTRL_LAB
  
  # ── Singleton filter ──────────────────────────────────────────────────────
  # Keep cells that are:
  #   (a) called Singleton by SNP demux
  #   (b) not a multi-donor barcode (no "+" in subject string)
  #   (c) assigned to a known subject (DONORS, SEVERE, or HC)
  is_singleton  <- tolower(as.character(md$SNP_call)) == "singleton"
  no_plus       <- !grepl("\\+", md$Subject_raw, perl = TRUE)
  ## PBMC CHANGE: whitelist now includes SEVERE subjects
  in_whitelist  <- md$Subject_fixed %in% c(DONORS, SEVERE, CTRL_LAB)
  
  # Also require the cell to have a valid HTO pool assignment
  ## PBMC CHANGE: exclude HTO doublets and negatives at the singleton stage
  valid_hto <- md$HTO_pool %in% c("Pool1", "Pool2")
  
  keep <- is_singleton & no_plus & in_whitelist & valid_hto
  
  log_status("FILTER",
             sprintf("Singleton + valid-HTO filter: %d/%d cells (%.1f%%)",
                     sum(keep), nrow(md), 100 * mean(keep)), "INFO")
  
  obj@meta.data <- md
  obj <- subset(obj, cells = rownames(obj@meta.data)[keep])
  obj <- join_if_multi(obj, stage_name = "JOIN_AFTER_SINGLETON_SUBSET")
  md  <- obj@meta.data
  
  # ── HTO-based TimePoint assignment ───────────────────────────────────────
  ## PBMC CHANGE: pool membership (HTO_pool) replaces lane-based logic
  ##   APECED patients: Pool1 = T1 = "Before", Pool2 = T2 = "After"
  ##   HC:              both pools → "Spike"  (spike-in control, no timepoint)
  ##   Severe:          Pool1 → "Severe_Pool1", Pool2 → "Severe_Pool2"
  
  assign_timepoint <- function(subject_fixed, hto_pool) {
    if (is.na(subject_fixed) || is.na(hto_pool)) return("Unknown")
    
    if (subject_fixed == CTRL_LAB)           return("Spike")
    
    if (subject_fixed %in% SEVERE) {
      if (hto_pool == "Pool1")               return("Severe_Pool1")
      if (hto_pool == "Pool2")               return("Severe_Pool2")
      return("Severe_Unknown")
    }
    
    if (subject_fixed %in% DONORS) {
      if (hto_pool == "Pool1")               return("Before")
      if (hto_pool == "Pool2")               return("After")
      return("Unknown")
    }
    
    return("Unknown")
  }
  
  md$TimePoint_fixed <- mapply(assign_timepoint,
                               md$Subject_fixed, md$HTO_pool,
                               USE.NAMES = FALSE)
  
  # ── Treatment label ───────────────────────────────────────────────────────
  ## PBMC CHANGE: Treatment column expanded to cover all four groups
  md$Treatment <- dplyr::case_when(
    md$TimePoint_fixed == "Before"        ~ "Untreated",
    md$TimePoint_fixed == "After"         ~ "Treated",
    md$TimePoint_fixed == "Spike"         ~ "HC",
    md$TimePoint_fixed == "Severe_Pool1"  ~ "Severe",
    md$TimePoint_fixed == "Severe_Pool2"  ~ "Severe",
    TRUE                                  ~ "Unknown"
  )
  
  obj@meta.data <- md
  
  # Summary
  log_status("FILTER",
             sprintf("Final cell count: %s", format(ncol(obj), big.mark = ",")),
             "SUCCESS")
  tp_tbl <- table(obj$TimePoint_fixed)
  for (tp in names(tp_tbl))
    log_status("FILTER",
               sprintf("  TimePoint %-18s : %s cells", tp,
                       format(tp_tbl[[tp]], big.mark = ",")), "INFO")
  
  return(obj)
}

################################################################################
# STAGE 4: DIMENSIONAL REDUCTION  (unchanged)
################################################################################

stage4_dimension_reduction <- function(obj) {
  section_header("STAGE 4: DIMENSIONAL REDUCTION")
  
  obj <- join_if_multi(obj, stage_name = "JOIN_BEFORE_DIMRED")
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject"))
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  
  log_status("DIMRED", "Normalizing...",              "PROGRESS")
  obj <- NormalizeData(obj, verbose = FALSE)
  
  log_status("DIMRED", "Finding variable features...", "PROGRESS")
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 3000, verbose = FALSE)
  
  log_status("DIMRED", "Scaling...",                  "PROGRESS")
  obj <- ScaleData(obj, verbose = FALSE)
  
  log_status("DIMRED", "Running PCA...",               "PROGRESS")
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  
  log_status("DIMRED", "Computing neighbors and clusters...", "PROGRESS")
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj,  resolution = 0.5, verbose = FALSE)
  
  log_status("DIMRED", "Running UMAP...",              "PROGRESS")
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  
  log_status("DIMRED",
             sprintf("Identified %d clusters",
                     length(unique(obj$seurat_clusters))), "SUCCESS")
  return(obj)
}

################################################################################
# STAGE 5: QC FILTERING  (unchanged)
################################################################################

stage5_qc_filtering <- function(obj) {
  section_header("STAGE 5: QC FILTERING")
  
  if (!"percent.mt" %in% colnames(obj@meta.data))
    obj$percent.mt <- obj$percent.mito
  
  if (!"percent.ribo" %in% colnames(obj@meta.data)) {
    ribo <- grep("^RPL|^RPS", rownames(obj), perl = TRUE, value = TRUE)
    if (length(ribo) > 0) {
      counts_mat  <- GetAssayData(obj, slot = "counts")
      ribo_counts <- Matrix::colSums(counts_mat[ribo, , drop = FALSE])
      total_counts <- Matrix::colSums(counts_mat)
      obj$percent.ribo <- (ribo_counts / pmax(total_counts, 1)) * 100
    }
  }
  
  pre_filter   <- ncol(obj)
  obj_filtered <- subset(
    obj,
    subset = nCount_RNA   > QC_THRESHOLDS$nCount_min   &
      nCount_RNA   < QC_THRESHOLDS$nCount_max   &
      nFeature_RNA > QC_THRESHOLDS$nFeature_min &
      nFeature_RNA < QC_THRESHOLDS$nFeature_max &
      percent.mt   < QC_THRESHOLDS$percent_mito_max &
      percent.ribo < QC_THRESHOLDS$percent_ribo_max &
      percent.hb   < QC_THRESHOLDS$percent_hb_max
  )
  
  obj_filtered <- join_if_multi(obj_filtered, stage_name = "JOIN_AFTER_QC")
  post_filter  <- ncol(obj_filtered)
  
  log_status("FILTER",
             sprintf("%s/%s cells retained (%.1f%%)",
                     format(post_filter, big.mark = ","),
                     format(pre_filter,  big.mark = ","),
                     100 * post_filter / pre_filter), "SUCCESS")
  return(obj_filtered)
}

################################################################################
# STAGE 6: CD45 GATING  (unchanged)
################################################################################

stage6_cd45_gating <- function(obj) {
  section_header("STAGE 6: CD45/EPCAM GATING")
  
  rn         <- rownames(obj)
  cd45_gene  <- if ("PTPRC" %in% rn) "PTPRC" else if ("CD45" %in% rn) "CD45" else NA_character_
  epcam_gene <- if ("EPCAM" %in% rn) "EPCAM" else NA_character_
  
  obj$CD45_expr  <- if (!is.na(cd45_gene))  FetchData(obj, cd45_gene)[, 1]  else 0
  obj$EPCAM_expr <- if (!is.na(epcam_gene)) FetchData(obj, epcam_gene)[, 1] else 0
  
  thr_kmeans <- function(x) {
    x0 <- x[x > 0]
    if (length(unique(x0)) < 2) return(0)
    km <- kmeans(x0, centers = 2, nstart = 10)
    m  <- tapply(x0, km$cluster, mean)
    mean(sort(m))
  }
  
  thr_cd45  <- thr_kmeans(obj$CD45_expr)
  thr_epcam <- thr_kmeans(obj$EPCAM_expr)
  obj$CD45pos <- (obj$CD45_expr > thr_cd45) & (obj$EPCAM_expr <= thr_epcam)
  
  log_status("GATE",
             sprintf("CD45+ cells: %s (%.1f%%)",
                     format(sum(obj$CD45pos), big.mark = ","),
                     100 * sum(obj$CD45pos) / ncol(obj)), "SUCCESS")
  return(obj)
}

################################################################################
# STAGE 7: SINGLER ANNOTATION  (unchanged)
################################################################################

stage7_singler_annotation <- function(obj) {
  section_header("STAGE 7: SINGLER REFERENCE ANNOTATION")
  
  if (!requireNamespace("SingleR",  quietly = TRUE)) {
    log_status("SINGLER", "SingleR not installed - skipping", "WARNING")
    return(obj)
  }
  if (!requireNamespace("celldex", quietly = TRUE)) {
    log_status("SINGLER", "celldex not installed - skipping", "WARNING")
    return(obj)
  }
  
  obj         <- check_and_join_layers(obj, "SINGLER")
  data_layers <- grep("^data", SeuratObject::Layers(obj[["RNA"]]), value = TRUE)
  
  if (length(data_layers) > 0) {
    counts_for_singler <- SeuratObject::LayerData(obj[["RNA"]],
                                                  layer = data_layers[1])
  } else {
    log_status("SINGLER", "No normalised data – normalising now", "INFO")
    obj         <- NormalizeData(obj, verbose = FALSE)
    data_layers <- grep("^data", SeuratObject::Layers(obj[["RNA"]]), value = TRUE)
    counts_for_singler <- SeuratObject::LayerData(obj[["RNA"]],
                                                  layer = data_layers[1])
  }
  
  log_status("SINGLER",
             sprintf("Matrix: %d genes × %d cells",
                     nrow(counts_for_singler), ncol(counts_for_singler)), "INFO")
  
  # Monaco Reference
  log_status("SINGLER", "Loading Monaco Immune reference...", "PROGRESS")
  tryCatch({
    ref_monaco       <- celldex::MonacoImmuneData()
    pred_monaco_main <- SingleR(test = counts_for_singler, ref = ref_monaco,
                                labels = ref_monaco$label.main,
                                BPPARAM = BiocParallel::SerialParam())
    pred_monaco_fine <- SingleR(test = counts_for_singler, ref = ref_monaco,
                                labels = ref_monaco$label.fine,
                                BPPARAM = BiocParallel::SerialParam())
    
    obj$Monaco_main         <- pred_monaco_main$labels
    obj$Monaco_main_pruned  <- pred_monaco_main$pruned.labels
    obj$Monaco_main_scores  <- pred_monaco_main$scores[
      cbind(seq_len(nrow(pred_monaco_main)),
            match(pred_monaco_main$labels, colnames(pred_monaco_main$scores)))]
    obj$Monaco_fine         <- pred_monaco_fine$labels
    obj$Monaco_fine_pruned  <- pred_monaco_fine$pruned.labels
    obj$Monaco_fine_scores  <- pred_monaco_fine$scores[
      cbind(seq_len(nrow(pred_monaco_fine)),
            match(pred_monaco_fine$labels, colnames(pred_monaco_fine$scores)))]
    
    log_status("SINGLER",
               sprintf("Monaco main: %d labels, %.1f%% retained after pruning",
                       length(unique(obj$Monaco_main)),
                       100 * sum(!is.na(obj$Monaco_main_pruned)) / ncol(obj)),
               "SUCCESS")
  }, error = function(e)
    log_status("SINGLER",
               sprintf("Monaco failed: %s", e$message), "WARNING"))
  
  # HPCA Reference
  log_status("SINGLER", "Loading HPCA reference...", "PROGRESS")
  tryCatch({
    ref_hpca       <- celldex::HumanPrimaryCellAtlasData()
    pred_hpca_main <- SingleR(test = counts_for_singler, ref = ref_hpca,
                              labels = ref_hpca$label.main,
                              BPPARAM = BiocParallel::SerialParam())
    pred_hpca_fine <- SingleR(test = counts_for_singler, ref = ref_hpca,
                              labels = ref_hpca$label.fine,
                              BPPARAM = BiocParallel::SerialParam())
    
    obj$HPCA_main         <- pred_hpca_main$labels
    obj$HPCA_main_pruned  <- pred_hpca_main$pruned.labels
    obj$HPCA_main_scores  <- pred_hpca_main$scores[
      cbind(seq_len(nrow(pred_hpca_main)),
            match(pred_hpca_main$labels, colnames(pred_hpca_main$scores)))]
    obj$HPCA_fine         <- pred_hpca_fine$labels
    obj$HPCA_fine_pruned  <- pred_hpca_fine$pruned.labels
    obj$HPCA_fine_scores  <- pred_hpca_fine$scores[
      cbind(seq_len(nrow(pred_hpca_fine)),
            match(pred_hpca_fine$labels, colnames(pred_hpca_fine$scores)))]
    
    log_status("SINGLER",
               sprintf("HPCA main: %d labels, %.1f%% retained after pruning",
                       length(unique(obj$HPCA_main)),
                       100 * sum(!is.na(obj$HPCA_main_pruned)) / ncol(obj)),
               "SUCCESS")
  }, error = function(e)
    log_status("SINGLER",
               sprintf("HPCA failed: %s", e$message), "WARNING"))
  
  singler_dir  <- file.path(OUT_DIR, "singler_results")
  dir.create(singler_dir, showWarnings = FALSE)
  singler_cols <- grep("^(Monaco|HPCA)_", colnames(obj@meta.data), value = TRUE)
  if (length(singler_cols) > 0) {
    singler_df         <- obj@meta.data[, singler_cols, drop = FALSE]
    singler_df$cell_id <- rownames(singler_df)
    fwrite(singler_df, file.path(singler_dir, "singler_per_cell_annotations.csv"))
    log_status("SINGLER", "SAVED: singler_per_cell_annotations.csv", "SUCCESS")
  }
  
  return(obj)
}

################################################################################
# STAGE 8: CLUSTER ANNOTATION  (unchanged)
################################################################################

stage8_cluster_annotation <- function(obj) {
  section_header("STAGE 8: CLUSTER ANNOTATION")
  
  ## PBMC CHANGE: EpithelialScore, lineage calling, and all tissue-specific
  ## cluster stats removed. PBMC contains only immune cells so these metrics
  ## are not meaningful. Stage 8 now only: finds marker genes (with caching),
  ## aggregates SingleR labels per cluster, builds the annotation table, and
  ## adds cluster_annotation to metadata.
  
  markers_dir       <- file.path(OUT_DIR, "markers")
  dir.create(markers_dir, showWarnings = FALSE, recursive = TRUE)
  
  markers_all_file  <- file.path(markers_dir, "FindAllMarkers_all_results.csv")
  markers_top5_file <- file.path(markers_dir, "FindAllMarkers_top5_per_cluster.csv")
  markers_sig_file  <- file.path(markers_dir, "FindAllMarkers_significant_only.csv")
  
  obj <- join_if_multi(obj, stage_name = "JOIN_BEFORE_CLUSTER_ANNOT")
  obj <- check_and_join_layers(obj, "ANNOT")
  
  obj$cluster <- factor(obj@meta.data[["seurat_clusters"]])
  md          <- obj@meta.data
  md$cluster  <- as.character(obj$cluster)
  
  # ── Cluster cell counts ───────────────────────────────────────────────────
  cluster_sizes <- md %>%
    group_by(cluster) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    mutate(pct_of_total = round(100 * n_cells / ncol(obj), 2))
  
  # ── Marker genes (with caching) ───────────────────────────────────────────
  load_cached <- !ANNOTATION_CONFIG$force_recompute &&
    file.exists(markers_all_file) &&
    file.exists(markers_top5_file)
  
  if (load_cached) {
    log_status("ANNOT", "Loading cached FindAllMarkers results...", "INFO")
    all_markers <- fread(markers_all_file)
    markers_df  <- fread(markers_top5_file)
    log_status("ANNOT",
               sprintf("Loaded %s markers from cache",
                       format(nrow(all_markers), big.mark = ",")), "SUCCESS")
  } else {
    log_status("ANNOT", "Finding marker genes (this may take several minutes)...",
               "PROGRESS")
    all_markers <- tryCatch(
      FindAllMarkers(object          = obj,
                     only.pos        = ANNOTATION_CONFIG$only_pos,
                     min.pct         = ANNOTATION_CONFIG$min_pct,
                     logfc.threshold = ANNOTATION_CONFIG$logfc_threshold,
                     return.thresh   = ANNOTATION_CONFIG$p_adj_max,
                     verbose         = FALSE),
      error = function(e) {
        log_status("ANNOT", sprintf("Marker finding error: %s", e$message), "WARNING")
        NULL
      }
    )
    
    if (!is.null(all_markers) && nrow(all_markers) > 0) {
      all_markers$cluster <- as.character(all_markers$cluster)
      fwrite(all_markers, markers_all_file)
      log_status("ANNOT",
                 sprintf("SAVED: %s (%s markers)",
                         basename(markers_all_file),
                         format(nrow(all_markers), big.mark = ",")), "SUCCESS")
      
      sig_markers <- all_markers %>% filter(p_val_adj <= ANNOTATION_CONFIG$p_adj_max)
      fwrite(sig_markers, markers_sig_file)
      
      markers_df <- all_markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        slice_head(n = 5) %>%
        summarise(top_markers = paste0(gene,             collapse = ", "),
                  top_logFC   = paste0(round(avg_log2FC, 2), collapse = ", "),
                  .groups = "drop")
      fwrite(markers_df, markers_top5_file)
      
      for (clust in unique(all_markers$cluster)) {
        fwrite(all_markers %>% filter(cluster == clust),
               file.path(markers_dir, sprintf("cluster_%s_markers.csv", clust)))
      }
      log_status("ANNOT",
                 sprintf("SAVED: per-cluster files for %d clusters",
                         length(unique(all_markers$cluster))), "SUCCESS")
    } else {
      markers_df <- data.frame(cluster     = character(),
                               top_markers = character(),
                               top_logFC   = character())
      log_status("ANNOT", "No significant markers found", "WARNING")
    }
  }
  
  # ── Aggregate SingleR labels per cluster (modal label) ───────────────────
  get_mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    if (length(ux) == 0) return(NA_character_)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  monaco_df <- if (any(c("Monaco_main", "Monaco_main_pruned") %in% colnames(md))) {
    md %>% group_by(cluster) %>%
      summarise(monaco_main        = get_mode(Monaco_main),
                monaco_main_pruned = get_mode(Monaco_main_pruned),
                .groups = "drop")
  } else {
    data.frame(cluster            = unique(md$cluster),
               monaco_main        = NA_character_,
               monaco_main_pruned = NA_character_)
  }
  
  hpca_df <- if (any(c("HPCA_main", "HPCA_main_pruned") %in% colnames(md))) {
    md %>% group_by(cluster) %>%
      summarise(hpca_main        = get_mode(HPCA_main),
                hpca_main_pruned = get_mode(HPCA_main_pruned),
                .groups = "drop")
  } else {
    data.frame(cluster          = unique(md$cluster),
               hpca_main        = NA_character_,
               hpca_main_pruned = NA_character_)
  }
  
  # ── Build and save annotation table ──────────────────────────────────────
  final_annot <- cluster_sizes %>%
    mutate(cluster = as.character(cluster)) %>%
    left_join(markers_df %>% mutate(cluster = as.character(cluster)), by = "cluster") %>%
    left_join(monaco_df  %>% mutate(cluster = as.character(cluster)), by = "cluster") %>%
    left_join(hpca_df    %>% mutate(cluster = as.character(cluster)), by = "cluster")
  
  final_annot$annotation_summary <- with(final_annot, paste0(
    "Monaco: ", ifelse(is.na(monaco_main), "NA", monaco_main),
    " | HPCA: ", ifelse(is.na(hpca_main),  "NA", hpca_main)
  ))
  
  fwrite(final_annot, file.path(OUT_DIR, "cluster_annotations_comprehensive.csv"))
  log_status("ANNOT", "SAVED: cluster_annotations_comprehensive.csv", "SUCCESS")
  
  annot_lookup <- setNames(final_annot$annotation_summary, final_annot$cluster)
  obj@meta.data$cluster_annotation <-
    annot_lookup[as.character(obj@meta.data$cluster)]
  
  log_status("ANNOT",
             sprintf("Added cluster_annotation: %d/%d cells",
                     sum(!is.na(obj@meta.data$cluster_annotation)), ncol(obj)),
             "SUCCESS")
  return(obj)
}

################################################################################
# STAGE 9: VISUALIZATION
## PBMC CHANGE:
# STAGE 9: VISUALIZATION
## PBMC CHANGE:
##   1. subject_pal expanded to include all 4 Severe subjects + HC
##   2. Added HTO pool UMAP plot
##   3. TimePoint comparison plot uses all six TimePoint_fixed levels
################################################################################

stage9_generate_plots <- function(obj) {
  section_header("STAGE 9: GENERATING VISUALIZATIONS")
  
  plots_dir        <- file.path(OUT_DIR, "plots")
  doublet_plots_dir <- file.path(plots_dir, "doublet_detection")
  dir.create(plots_dir,         showWarnings = FALSE)
  dir.create(doublet_plots_dir, showWarnings = FALSE)
  
  ## PBMC CHANGE: palette now covers APECED patients, all 4 Severe subjects, and HC
  subject_pal <- c(
    # APECED patients
    "P0001171" = "#F06565",
    "P0001352" = "#9ACD32",
    "P0001358" = "#00B5D8",
    "P0002345" = "#B085F5",
    # Severe patients
    "P0001140" = "#FF8C00",
    "P0011505" = "#FF4500",
    "P0005819" = "#DC143C",
    "P0014164" = "#8B0000",
    # Healthy controls (collapsed label)
    "HC"       = "#333333"
  )
  
  ## PBMC CHANGE: assign cluster alias here since Stage 9 runs before Stage 8.
  ## seurat_clusters is always present after FindClusters(); this just creates
  ## the factor alias used throughout the plot calls below.
  if (!"cluster" %in% colnames(obj@meta.data))
    obj$cluster <- factor(obj@meta.data[["seurat_clusters"]])
  
  log_status("PLOT", "Generating UMAP plots...", "PROGRESS")
  
  # Subject
  p1 <- DimPlot(obj, reduction = "umap", group.by = "Subject_fixed",
                cols = subject_pal, pt.size = 0.2, raster = TRUE) +
    ggtitle("UMAP coloured by Subject")
  save_plot_quietly(file.path(plots_dir, "UMAP_by_Subject.png"), p1,
                    width = 9, height = 6.5, dpi = 300)
  
  # TimePoint
  if ("TimePoint_fixed" %in% colnames(obj@meta.data)) {
    ## PBMC CHANGE: palette covers all six TimePoint levels
    tp_pal <- c("Before"       = "#4575B4",
                "After"        = "#D73027",
                "Spike"        = "#333333",
                "Severe_Pool1" = "#FC8D59",
                "Severe_Pool2" = "#FEE090",
                "Unknown"      = "#AAAAAA")
    p2 <- DimPlot(obj, reduction = "umap", group.by = "TimePoint_fixed",
                  cols = tp_pal, pt.size = 0.2, raster = TRUE) +
      ggtitle("UMAP coloured by TimePoint")
    save_plot_quietly(file.path(plots_dir, "UMAP_by_TimePoint.png"), p2,
                      width = 8, height = 6.5, dpi = 300)
  }
  
  # Cluster
  p3 <- DimPlot(obj, reduction = "umap", group.by = "cluster",
                pt.size = 0.2, raster = TRUE, label = TRUE) +
    ggtitle("UMAP coloured by Cluster")
  save_plot_quietly(file.path(plots_dir, "UMAP_by_Cluster.png"), p3,
                    width = 9, height = 7, dpi = 300)
  
  # Subject split by Lane (only if Lane column was populated by SNP demux)
  if ("Lane" %in% colnames(obj@meta.data)) {
    p4 <- DimPlot(obj, reduction = "umap", group.by = "Subject_fixed",
                  split.by = "Lane", cols = subject_pal, ncol = 3,
                  pt.size = 0.2, raster = TRUE) +
      ggtitle("UMAP by Subject, split by Lane")
    save_plot_quietly(file.path(plots_dir, "UMAP_by_Subject_split_Lane.png"), p4,
                      width = 16, height = 10, dpi = 300)
  } else {
    log_status("PLOT", "Skipping Lane split UMAP – Lane column not present", "WARNING")
  }
  
  ## PBMC CHANGE: HTO pool UMAP
  if ("HTO_pool" %in% colnames(obj@meta.data)) {
    pool_pal <- c("Pool1"       = "#1F78B4",
                  "Pool2"       = "#E31A1C",
                  "HTO_Doublet" = "#6A3D9A",
                  "HTO_Negative"= "#B15928",
                  "HTO_Unknown" = "#AAAAAA")
    p_hto <- DimPlot(obj, reduction = "umap", group.by = "HTO_pool",
                     cols = pool_pal, pt.size = 0.2, raster = TRUE) +
      ggtitle("UMAP coloured by HTO Pool (HTODemux)")
    save_plot_quietly(file.path(plots_dir, "UMAP_by_HTO_pool.png"), p_hto,
                      width = 8, height = 6.5, dpi = 300)
  }
  
  # ── Doublet detection plots (unchanged logic) ─────────────────────────────
  log_status("PLOT", "Generating doublet detection plots...", "PROGRESS")
  
  doublet_colors <- c("singlet" = "#2E7D32", "doublet" = "#D32F2F",
                      "unknown" = "#757575")
  
  if ("doublet_consensus" %in% colnames(obj@meta.data)) {
    p_consensus <- DimPlot(obj, reduction = "umap", group.by = "doublet_consensus",
                           cols = doublet_colors, pt.size = 1, raster = TRUE) +
      ggtitle("UMAP: Consensus Doublet Detection")
    save_plot_quietly(file.path(doublet_plots_dir, "UMAP_doublet_consensus.png"),
                      p_consensus, width = 9, height = 7, dpi = 300)
  }
  
  doublet_tools <- c("DoubletDetection_Droplet", "scDblFinder_Droplet",
                     "scds_Droplet",             "scrublet_Droplet")
  tool_names    <- c("DoubletDetection", "scDblFinder", "scds", "scrublet")
  available_tools <- doublet_tools[doublet_tools %in% colnames(obj@meta.data)]
  
  if (length(available_tools) > 0) {
    plot_list <- list()
    for (i in seq_along(available_tools)) {
      tool_col  <- available_tools[i]
      tool_name <- tool_names[i]
      p_tool <- DimPlot(obj, reduction = "umap", group.by = tool_col,
                        cols = doublet_colors, pt.size = 1, raster = TRUE) +
        ggtitle(tool_name) +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.5, face = "bold"))
      plot_list[[tool_name]] <- p_tool
      save_plot_quietly(
        file.path(doublet_plots_dir, sprintf("UMAP_doublet_%s.png", tool_name)),
        p_tool, width = 7, height = 6, dpi = 300)
    }
    if (length(plot_list) == 4) {
      p_combined <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(title = "Doublet Detection: Comparison Across Tools",
                        theme = theme(plot.title =
                                        element_text(size = 16, face = "bold",
                                                     hjust = 0.5)))
      save_plot_quietly(
        file.path(doublet_plots_dir, "UMAP_doublet_all_tools_comparison.png"),
        p_combined, width = 14, height = 12, dpi = 300)
    }
  }
  
  score_cols       <- c("scDblFinder_Score", "scds_Score", "scrublet_Score")
  available_scores <- score_cols[score_cols %in% colnames(obj@meta.data)]
  if (length(available_scores) > 0) {
    score_plots <- list()
    for (score_col in available_scores) {
      tool_name <- sub("_Score", "", score_col)
      score_plots[[tool_name]] <- FeaturePlot(
        obj, features = score_col, reduction = "umap",
        pt.size = 1, raster = TRUE) +
        scale_color_viridis_c(option = "magma", name = "Score") +
        ggtitle(paste(tool_name, "Score")) +
        theme(legend.position = "right",
              plot.title = element_text(hjust = 0.5, face = "bold"))
    }
    save_plot_quietly(
      file.path(doublet_plots_dir, "UMAP_doublet_scores_comparison.png"),
      wrap_plots(score_plots, ncol = 2) +
        plot_annotation(title = "Doublet Detection Scores by Tool",
                        theme = theme(plot.title =
                                        element_text(size = 16, face = "bold",
                                                     hjust = 0.5))),
      width = 14, height = 10, dpi = 300)
  }
  
  if ("doublet_consensus" %in% colnames(obj@meta.data)) {
    doublet_summary <- as.data.frame(obj@meta.data) %>%
      select(contains("Droplet")) %>%
      pivot_longer(cols = everything(), names_to = "Tool", values_to = "Call") %>%
      filter(!is.na(Call)) %>%
      group_by(Tool, Call) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Tool) %>%
      mutate(Percentage = 100 * Count / sum(Count),
             Tool = sub("_Droplet", "", Tool))
    
    p_summary <- ggplot(doublet_summary,
                        aes(x = Tool, y = Percentage, fill = Call)) +
      geom_col(position = "stack") +
      geom_text(aes(label = sprintf("%.1f%%", Percentage)),
                position = position_stack(vjust = 0.5),
                size = 3.5, color = "white") +
      scale_fill_manual(values = doublet_colors) +
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Doublet Detection Summary Across Tools",
           x = "Tool", y = "Percentage of Cells", fill = "Classification")
    save_plot_quietly(
      file.path(doublet_plots_dir, "Barplot_doublet_summary.png"),
      p_summary, width = 10, height = 6, dpi = 300)
  }
  
  # ── QC feature plots ──────────────────────────────────────────────────────
  log_status("PLOT", "Generating QC plots...", "PROGRESS")
  
  p6 <- FeaturePlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                    reduction = "umap", order = TRUE, raster = TRUE, ncol = 3) &
    NoLegend()
  save_plot_quietly(file.path(plots_dir, "UMAP_QC_metrics.png"), p6,
                    width = 15, height = 5, dpi = 300)
  
  markers <- intersect(c("EPCAM","KRT8","KRT18","PTPRC","CD3D","CD14"),
                       rownames(obj))
  if (length(markers) > 0) {
    p7 <- FeaturePlot(obj, features = markers, reduction = "umap",
                      order = TRUE, raster = TRUE, ncol = 3,
                      min.cutoff = "q5", max.cutoff = "q95") & NoLegend()
    save_plot_quietly(file.path(plots_dir, "UMAP_key_markers.png"), p7,
                      width = 15, height = 10, dpi = 300)
  }
  
  # SingleR plots
  for (col in c("Monaco_main","Monaco_main_pruned","HPCA_main","HPCA_main_pruned")) {
    if (col %in% colnames(obj@meta.data)) {
      p_ref <- DimPlot(obj, reduction = "umap", group.by = col,
                       pt.size = 0.2, raster = TRUE,
                       label = TRUE, repel = TRUE) +
        ggtitle(paste("UMAP:", col))
      save_plot_quietly(file.path(plots_dir, paste0("UMAP_", col, ".png")),
                        p_ref, width = 10, height = 8, dpi = 300)
    }
  }
  
  # Violin QC
  v1 <- VlnPlot(obj,
                features  = c("nFeature_RNA","nCount_RNA","percent.mt"),
                group.by  = "Subject_fixed", pt.size = 0, ncol = 3)
  save_plot_quietly(file.path(plots_dir, "Violin_QC_by_Subject.png"), v1,
                    width = 16, height = 5, dpi = 300)
  
  # Cluster composition Before vs After (APECED patients only)
  if ("TimePoint_fixed" %in% colnames(obj@meta.data)) {
    comp_data <- as.data.frame(obj@meta.data) %>%
      mutate(cluster = as.character(cluster)) %>%
      filter(TimePoint_fixed %in% c("Before", "After")) %>%
      group_by(TimePoint_fixed, cluster) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(TimePoint_fixed) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
    
    p9 <- ggplot(comp_data,
                 aes(x = cluster, y = frac, fill = TimePoint_fixed)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7) +
      scale_y_continuous(labels = percent_format()) +
      theme_bw(12) +
      labs(title = "Cluster composition: Before vs After (APECED patients)",
           x = "Cluster", y = "Fraction", fill = "TimePoint")
    save_plot_quietly(
      file.path(plots_dir, "Barplot_cluster_composition_BeforeAfter.png"),
      p9, width = 12, height = 6, dpi = 300)
    
    ## PBMC CHANGE: additional plot comparing all four treatment groups
    comp_all <- as.data.frame(obj@meta.data) %>%
      mutate(cluster = as.character(cluster)) %>%
      filter(TimePoint_fixed %in% c("Before","After",
                                    "Severe_Pool1","Severe_Pool2","Spike")) %>%
      group_by(TimePoint_fixed, cluster) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(TimePoint_fixed) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
    
    p9b <- ggplot(comp_all,
                  aes(x = cluster, y = frac, fill = TimePoint_fixed)) +
      geom_col(position = position_dodge(width = 0.9), width = 0.8) +
      scale_y_continuous(labels = percent_format()) +
      scale_fill_manual(values = c("Before"       = "#4575B4",
                                   "After"        = "#D73027",
                                   "Spike"        = "#333333",
                                   "Severe_Pool1" = "#FC8D59",
                                   "Severe_Pool2" = "#FEE090")) +
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Cluster composition: all groups",
           x = "Cluster", y = "Fraction", fill = "Group")
    save_plot_quietly(
      file.path(plots_dir, "Barplot_cluster_composition_AllGroups.png"),
      p9b, width = 14, height = 6, dpi = 300)
  }
  
  # Cell count heatmap Lane × Subject (only if Lane column present)
  if ("Lane" %in% colnames(obj@meta.data)) {
    count_data           <- as.data.frame(table(obj$Lane, obj$Subject_fixed))
    colnames(count_data) <- c("Lane", "Subject", "Count")
    p10 <- ggplot(count_data, aes(x = Subject, y = Lane, fill = Count)) +
      geom_tile() +
      geom_text(aes(label = scales::comma(Count)), size = 3) +
      scale_fill_gradient(low = "white", high = "#3B82F6") +
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Cell counts (Lane × Subject)", x = "Subject", y = "Lane")
    save_plot_quietly(file.path(plots_dir, "Heatmap_counts_Lane_vs_Subject.png"),
                      p10, width = 10, height = 5.5, dpi = 300)
  } else {
    log_status("PLOT", "Skipping Lane heatmap – Lane column not present", "WARNING")
  }
  
  # ============================================================
  # ANNOTATION BARPLOTS
  # ============================================================
  log_status("PLOT", "Generating annotation barplots...", "PROGRESS")
  
  annot_qc_dir <- file.path(plots_dir, "annotation_barplots")
  dir.create(annot_qc_dir, showWarnings = FALSE)
  
  annot_cols <- c("Monaco_main", "HPCA_main", "cluster_annotation")
  avail_annot <- annot_cols[annot_cols %in% colnames(obj@meta.data)]
  
  group_vars <- list(
    TimePoint = "TimePoint_fixed",
    Subject   = "Subject_fixed",
    Pool      = "HTO_pool"
  )
  
  for (annot_col in avail_annot) {
    annot_label <- sub("_main$|_annotation$", "", annot_col)
    
    for (gname in names(group_vars)) {
      gvar <- group_vars[[gname]]
      if (!gvar %in% colnames(obj@meta.data)) next
      
      df <- as.data.frame(obj@meta.data) %>%
        filter(!is.na(.data[[annot_col]]), !is.na(.data[[gvar]])) %>%
        group_by(.data[[gvar]], .data[[annot_col]]) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(.data[[gvar]]) %>%
        mutate(frac = n / sum(n)) %>%
        ungroup()
      colnames(df)[1:2] <- c("Group", "CellType")
      
      n_types <- length(unique(df$CellType))
      pal <- if (n_types <= 8)  RColorBrewer::brewer.pal(max(3, n_types), "Set2")[seq_len(n_types)] else
        if (n_types <= 12) RColorBrewer::brewer.pal(12, "Paired")[seq_len(n_types)] else
          scales::hue_pal()(n_types)
      
      p_bar <- ggplot(df, aes(x = Group, y = frac, fill = CellType)) +
        geom_col(position = "stack", width = 0.75) +
        scale_y_continuous(labels = percent_format(accuracy = 1)) +
        scale_fill_manual(values = setNames(pal, unique(df$CellType))) +
        theme_bw(12) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1),
              legend.title  = element_text(size = 10),
              legend.text   = element_text(size = 8),
              legend.key.size = unit(0.4, "cm")) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(title = sprintf("%s composition by %s", annot_label, gname),
             x = gname, y = "Fraction", fill = annot_label)
      
      fname <- sprintf("Barplot_%s_by_%s.png", annot_label, gname)
      save_plot_quietly(file.path(annot_qc_dir, fname), p_bar,
                        width = max(8, length(unique(df$Group)) * 1.2 + 4),
                        height = 7, dpi = 300)
    }
    
    # Per-subject stacked (one bar per subject, grouped by TimePoint colour border)
    if ("Subject_fixed" %in% colnames(obj@meta.data) &&
        "TimePoint_fixed" %in% colnames(obj@meta.data)) {
      df_subj <- as.data.frame(obj@meta.data) %>%
        filter(!is.na(.data[[annot_col]])) %>%
        group_by(Subject_fixed, TimePoint_fixed, .data[[annot_col]]) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(Subject_fixed) %>%
        mutate(frac = n / sum(n),
               subj_label = paste0(Subject_fixed, "
(", TimePoint_fixed, ")")) %>%
        ungroup()
      colnames(df_subj)[3] <- "CellType"
      
      n_types2 <- length(unique(df_subj$CellType))
      pal2 <- if (n_types2 <= 8)  RColorBrewer::brewer.pal(max(3, n_types2), "Set2")[seq_len(n_types2)] else
        if (n_types2 <= 12) RColorBrewer::brewer.pal(12, "Paired")[seq_len(n_types2)] else
          scales::hue_pal()(n_types2)
      
      p_subj <- ggplot(df_subj,
                       aes(x = subj_label, y = frac, fill = CellType)) +
        geom_col(position = "stack", width = 0.75) +
        scale_y_continuous(labels = percent_format(accuracy = 1)) +
        scale_fill_manual(values = setNames(pal2, unique(df_subj$CellType))) +
        theme_bw(12) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
              legend.text   = element_text(size = 8),
              legend.key.size = unit(0.4, "cm")) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(title = sprintf("%s composition per subject", annot_label),
             x = "Subject (TimePoint)", y = "Fraction", fill = annot_label)
      
      save_plot_quietly(
        file.path(annot_qc_dir, sprintf("Barplot_%s_perSubject.png", annot_label)),
        p_subj, width = max(10, length(unique(df_subj$subj_label)) * 0.9 + 4),
        height = 7, dpi = 300)
    }
  }
  
  # ============================================================
  # BATCH / QC PLOTS  (pool, lane, subject, timepoint overlays)
  # ============================================================
  log_status("PLOT", "Generating batch QC plots...", "PROGRESS")
  
  batch_dir <- file.path(plots_dir, "batch_qc")
  dir.create(batch_dir, showWarnings = FALSE)
  
  qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mito",
                  "percent.ribo", "percent.hb",
                  "ambient_rna_score", "library_complexity")
  qc_present <- qc_metrics[qc_metrics %in% colnames(obj@meta.data)]
  
  # ── 1. UMAP split by Pool (side-by-side panels) ───────────────────────────
  if ("HTO_pool" %in% colnames(obj@meta.data)) {
    pool_pal2 <- c("Pool1" = "#1F78B4", "Pool2" = "#E31A1C")
    p_pool_split <- DimPlot(obj, reduction = "umap", split.by = "HTO_pool",
                            group.by = "HTO_pool", cols = pool_pal2,
                            pt.size = 0.3, raster = TRUE) +
      ggtitle("UMAP split by HTO Pool")
    save_plot_quietly(file.path(batch_dir, "UMAP_Pool_split.png"),
                      p_pool_split, width = 12, height = 6, dpi = 300)
    
    # Pool overlaid on cluster
    p_pool_cluster <- DimPlot(obj, reduction = "umap", group.by = "HTO_pool",
                              cols = pool_pal2, pt.size = 0.3,
                              raster = TRUE, shuffle = TRUE) +
      ggtitle("UMAP: Pool1 vs Pool2 overlaid")
    save_plot_quietly(file.path(batch_dir, "UMAP_Pool_overlay.png"),
                      p_pool_cluster, width = 8, height = 6.5, dpi = 300)
    
    # Pool split by cluster – what fraction of each cluster is each pool?
    if ("cluster" %in% colnames(obj@meta.data)) {
      pool_clust <- as.data.frame(obj@meta.data) %>%
        filter(HTO_pool %in% c("Pool1","Pool2")) %>%
        group_by(cluster, HTO_pool) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(cluster) %>%
        mutate(frac = n / sum(n)) %>%
        ungroup()
      
      p_pool_clust_bar <- ggplot(pool_clust,
                                 aes(x = cluster, y = frac, fill = HTO_pool)) +
        geom_col(position = "stack", width = 0.75) +
        scale_fill_manual(values = pool_pal2) +
        scale_y_continuous(labels = percent_format(accuracy = 1)) +
        theme_bw(12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Pool composition per cluster",
             x = "Cluster", y = "Fraction", fill = "Pool")
      save_plot_quietly(file.path(batch_dir, "Barplot_Pool_per_cluster.png"),
                        p_pool_clust_bar, width = 10, height = 5, dpi = 300)
    }
  }
  
  # ── 2. UMAP split by Lane ─────────────────────────────────────────────────
  if ("Lane" %in% colnames(obj@meta.data)) {
    p_lane_split <- DimPlot(obj, reduction = "umap", split.by = "Lane",
                            group.by = "cluster", ncol = 3,
                            pt.size = 0.3, raster = TRUE, label = TRUE) +
      ggtitle("UMAP split by Lane") + NoLegend()
    save_plot_quietly(file.path(batch_dir, "UMAP_Lane_split.png"),
                      p_lane_split, width = 15, height = 10, dpi = 300)
    
    # Lane × Pool heatmap of cell counts
    if ("HTO_pool" %in% colnames(obj@meta.data)) {
      lane_pool <- as.data.frame(table(obj$Lane, obj$HTO_pool))
      colnames(lane_pool) <- c("Lane", "Pool", "Count")
      p_lp <- ggplot(lane_pool, aes(x = Pool, y = Lane, fill = Count)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = scales::comma(Count)), size = 3.5) +
        scale_fill_gradient(low = "white", high = "#3B82F6") +
        theme_bw(12) +
        labs(title = "Cell counts: Lane × Pool")
      save_plot_quietly(file.path(batch_dir, "Heatmap_Lane_x_Pool.png"),
                        p_lp, width = 6, height = 6, dpi = 300)
    }
  }
  
  # ── 3. QC violin by Pool ──────────────────────────────────────────────────
  if ("HTO_pool" %in% colnames(obj@meta.data) && length(qc_present) > 0) {
    v_pool <- VlnPlot(obj, features = qc_present, group.by = "HTO_pool",
                      pt.size = 0, ncol = min(4, length(qc_present)),
                      cols = c("Pool1" = "#1F78B4", "Pool2" = "#E31A1C")) +
      plot_annotation(title = "QC metrics by Pool")
    save_plot_quietly(file.path(batch_dir, "Violin_QC_by_Pool.png"),
                      v_pool,
                      width = min(4, length(qc_present)) * 3.5, height = 5,
                      dpi = 300)
  }
  
  # ── 4. QC violin by Lane ──────────────────────────────────────────────────
  if ("Lane" %in% colnames(obj@meta.data) && length(qc_present) > 0) {
    v_lane <- VlnPlot(obj, features = qc_present, group.by = "Lane",
                      pt.size = 0, ncol = min(4, length(qc_present))) +
      plot_annotation(title = "QC metrics by Lane")
    save_plot_quietly(file.path(batch_dir, "Violin_QC_by_Lane.png"),
                      v_lane,
                      width = min(4, length(qc_present)) * 3.5, height = 5,
                      dpi = 300)
  }
  
  # ── 5. QC feature plots on UMAP ──────────────────────────────────────────
  if (length(qc_present) > 0) {
    p_qc_feat <- FeaturePlot(obj, features = qc_present, reduction = "umap",
                             order = TRUE, raster = TRUE,
                             ncol = min(4, length(qc_present))) & NoLegend()
    save_plot_quietly(file.path(batch_dir, "UMAP_QC_all_metrics.png"),
                      p_qc_feat,
                      width = min(4, length(qc_present)) * 3.5,
                      height = ceiling(length(qc_present) / 4) * 3.5,
                      dpi = 300)
  }
  
  # ── 6. QC scatter: nCount vs nFeature coloured by Pool / mito ────────────
  md_qc <- as.data.frame(obj@meta.data)
  if (all(c("nCount_RNA","nFeature_RNA","percent.mito") %in% colnames(md_qc))) {
    
    p_scatter_mito <- ggplot(md_qc %>% slice_sample(n = min(50000, nrow(md_qc))),
                             aes(x = nCount_RNA, y = nFeature_RNA,
                                 colour = percent.mito)) +
      geom_point(size = 0.3, alpha = 0.5) +
      scale_colour_viridis_c(option = "magma", name = "% mito") +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      theme_bw(12) +
      labs(title = "nCount vs nFeature (colour = % mito)")
    save_plot_quietly(file.path(batch_dir, "Scatter_nCount_nFeature_mito.png"),
                      p_scatter_mito, width = 7, height = 5.5, dpi = 300)
    
    if ("HTO_pool" %in% colnames(md_qc)) {
      p_scatter_pool <- ggplot(
        md_qc %>% filter(HTO_pool %in% c("Pool1","Pool2")) %>%
          slice_sample(n = min(50000, nrow(.))),
        aes(x = nCount_RNA, y = nFeature_RNA, colour = HTO_pool)) +
        geom_point(size = 0.3, alpha = 0.4) +
        scale_colour_manual(values = c("Pool1"="#1F78B4","Pool2"="#E31A1C")) +
        scale_x_log10(labels = scales::comma) +
        scale_y_log10(labels = scales::comma) +
        theme_bw(12) +
        labs(title = "nCount vs nFeature by Pool", colour = "Pool")
      save_plot_quietly(file.path(batch_dir, "Scatter_nCount_nFeature_Pool.png"),
                        p_scatter_pool, width = 7, height = 5.5, dpi = 300)
    }
  }
  
  # ── 7. Cell counts per subject × pool (heatmap) ───────────────────────────
  if (all(c("Subject_fixed","HTO_pool") %in% colnames(obj@meta.data))) {
    subj_pool <- as.data.frame(table(obj$Subject_fixed, obj$HTO_pool))
    colnames(subj_pool) <- c("Subject", "Pool", "Count")
    p_sp <- ggplot(subj_pool, aes(x = Pool, y = Subject, fill = Count)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = scales::comma(Count)), size = 3.5) +
      scale_fill_gradient(low = "white", high = "#3B82F6") +
      theme_bw(12) +
      labs(title = "Cell counts: Subject × Pool")
    save_plot_quietly(file.path(batch_dir, "Heatmap_Subject_x_Pool.png"),
                      p_sp, width = 6, height = 8, dpi = 300)
  }
  
  # ── 8. Proportion of each subject per cluster (batch check) ───────────────
  if (all(c("Subject_fixed","cluster") %in% colnames(obj@meta.data))) {
    subj_clust <- as.data.frame(obj@meta.data) %>%
      group_by(cluster, Subject_fixed) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(cluster) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
    
    p_subj_clust <- ggplot(subj_clust,
                           aes(x = cluster, y = frac, fill = Subject_fixed)) +
      geom_col(position = "stack", width = 0.75) +
      scale_fill_manual(values = subject_pal) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Subject proportion per cluster (batch check)",
           x = "Cluster", y = "Fraction", fill = "Subject")
    save_plot_quietly(file.path(batch_dir, "Barplot_Subject_per_cluster.png"),
                      p_subj_clust, width = 12, height = 5, dpi = 300)
  }
  
  # ── 9. TimePoint proportion per cluster ───────────────────────────────────
  if (all(c("TimePoint_fixed","cluster") %in% colnames(obj@meta.data))) {
    tp_clust <- as.data.frame(obj@meta.data) %>%
      filter(!is.na(TimePoint_fixed)) %>%
      group_by(cluster, TimePoint_fixed) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(cluster) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
    
    tp_pal2 <- c("Before"       = "#4575B4", "After"        = "#D73027",
                 "Spike"        = "#333333", "Severe_Pool1" = "#FC8D59",
                 "Severe_Pool2" = "#FEE090", "Unknown"      = "#AAAAAA")
    
    p_tp_clust <- ggplot(tp_clust,
                         aes(x = cluster, y = frac, fill = TimePoint_fixed)) +
      geom_col(position = "stack", width = 0.75) +
      scale_fill_manual(values = tp_pal2) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "TimePoint proportion per cluster (batch check)",
           x = "Cluster", y = "Fraction", fill = "TimePoint")
    save_plot_quietly(file.path(batch_dir, "Barplot_TimePoint_per_cluster.png"),
                      p_tp_clust, width = 12, height = 5, dpi = 300)
  }
  
  log_status("PLOT", sprintf("Batch QC plots saved to: %s/", basename(batch_dir)),
             "SUCCESS")
  log_status("PLOT", sprintf("Annotation barplots saved to: %s/",
                             basename(annot_qc_dir)), "SUCCESS")
  
  log_status("PLOT", "All plots generated successfully", "SUCCESS")
  return(obj)
}

################################################################################
# MAIN PIPELINE
################################################################################

main <- function() {
  section_header("SINGLE-CELL RNA-SEQ ANALYSIS PIPELINE – PBMC (301-4)")
  
  start_time <- Sys.time()
  log_status("INIT", sprintf("Project: %s", PROJECT),      "INFO")
  log_status("INIT", sprintf("Output:  %s", OUT_DIR),      "INFO")
  log_status("INIT", sprintf("APECED donors:  %s", paste(DONORS, collapse = ", ")), "INFO")
  log_status("INIT", sprintf("Severe patients: %s", paste(SEVERE, collapse = ", ")), "INFO")
  log_status("INIT", sprintf("HTO Pool1: %s  |  Pool2: %s", HTO_POOL1, HTO_POOL2), "INFO")
  
  obj_raw <- stage1_load_and_integrate()
  obj_raw <- stage2_add_metadata(obj_raw)
  
  raw_file <- file.path(OUT_DIR, "object_raw_with_metadata.rds")
  log_status("SAVE", "Saving raw object with metadata...", "PROGRESS")
  saveRDS(obj_raw, raw_file)
  log_status("SAVE",
             sprintf("Saved: %s (%.1f MB)", basename(raw_file),
                     file.size(raw_file) / 1024^2), "SUCCESS")
  
  obj_clean    <- stage3_filter_and_harmonize(obj_raw)
  obj_clean    <- stage4_dimension_reduction(obj_clean)
  obj_filtered <- stage5_qc_filtering(obj_clean)
  
  log_status("DIMRED", "Re-running dimension reduction on filtered cells...", "PROGRESS")
  obj_filtered <- stage4_dimension_reduction(obj_filtered)
  
  ## PBMC CHANGE: Stage 6 CD45 gating removed – not applicable for PBMC
  obj_filtered <- stage7_singler_annotation(obj_filtered)
  ## PBMC CHANGE: plots run before marker finding so UMAPs are available
  ##   immediately; FindAllMarkers is the slow step and runs last.
  obj_filtered <- stage9_generate_plots(obj_filtered)
  obj_filtered <- stage8_cluster_annotation(obj_filtered)
  
  final_file <- file.path(OUT_DIR, "object_final_annotated.rds")
  log_status("SAVE", "Saving final annotated object...", "PROGRESS")
  saveRDS(obj_filtered, final_file)
  log_status("SAVE",
             sprintf("Saved: %s (%.1f MB)", basename(final_file),
                     file.size(final_file) / 1024^2), "SUCCESS")
  
  # Doublet summary CSV
  if ("doublet_consensus" %in% colnames(obj_filtered@meta.data)) {
    doublet_summary_df <- obj_filtered@meta.data %>%
      select(contains("Droplet"), contains("_Score"),
             doublet_consensus, doublet_confidence,
             n_doublet_calls, n_tools_called) %>%
      mutate(cell_id = rownames(obj_filtered@meta.data))
    fwrite(doublet_summary_df,
           file.path(OUT_DIR, "doublet_detection_summary.csv"))
  }
  
  section_header("PIPELINE SUMMARY")
  
  end_time <- Sys.time()
  elapsed  <- difftime(end_time, start_time, units = "mins")
  
  cat(sprintf("📊 ANALYSIS COMPLETE\n"))
  cat(sprintf("   Total runtime: %.1f minutes\n", as.numeric(elapsed)))
  cat(sprintf("\n📁 OUTPUT FILES:\n"))
  cat(sprintf("   1. Raw object:   %s\n",  basename(raw_file)))
  cat(sprintf("   2. Final object: %s\n",  basename(final_file)))
  cat(sprintf("   3. Annotations:  cluster_annotations_comprehensive.csv\n"))
  cat(sprintf("   4. Markers:      markers/ directory\n"))
  cat(sprintf("   5. SingleR:      singler_results/\n"))
  cat(sprintf("   6. Plots:        plots/ directory\n"))
  cat(sprintf("\n📈 DATA SUMMARY:\n"))
  cat(sprintf("   Raw cells:      %s\n", format(ncol(obj_raw),      big.mark = ",")))
  cat(sprintf("   Filtered cells: %s\n", format(ncol(obj_filtered), big.mark = ",")))
  cat(sprintf("   Clusters:       %d\n", length(unique(obj_filtered$cluster))))
  cat(sprintf("   Genes:          %s\n", format(nrow(obj_filtered), big.mark = ",")))
  
  
  if ("TimePoint_fixed" %in% colnames(obj_filtered@meta.data)) {
    cat(sprintf("\n⏱️  GROUP DISTRIBUTION:\n"))
    tp_table <- table(obj_filtered$TimePoint_fixed)
    for (tp in names(tp_table))
      cat(sprintf("   %-18s : %s cells (%.1f%%)\n",
                  tp,
                  format(tp_table[tp], big.mark = ","),
                  100 * tp_table[tp] / ncol(obj_filtered)))
  }
  
  cat(sprintf("\n✅ All outputs saved to: %s\n", OUT_DIR))
  cat(sprintf("💡 TIP: Set ANNOTATION_CONFIG$force_recompute = TRUE to regenerate markers\n\n"))
  
  log_status("COMPLETE", "Pipeline finished successfully!", "SUCCESS")
  invisible(obj_filtered)
}

################################################################################
# RUN
################################################################################

if (!interactive()) {
  main()
}
