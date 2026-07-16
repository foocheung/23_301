#!/usr/bin/env Rscript
# =============================================================================
# Cascade recovery analysis: PBMC disease phenotype -> PBMC treatment recovery
#                              -> Duodenum treatment recovery
#
# Answers three questions using the same file layout / naming conventions as
# "qmd_step3_pathway_strict_then_directional":
#
#   Q1) PBMC HC vs Mild (pooled): which pathways are FDR-significant
#       ("disease phenotype"), split Up vs Down, per Monaco fine cell type.
#
#   Q2) PBMC MildTreated vs MildUntreated (paired): of the Q1 disease-phenotype
#       pathways, which move back toward HC ("recovered in PBMC") vs not,
#       based on DIRECTION ONLY (no significance filter on the treatment stat).
#
#   Q3) Duodenum Treated vs Untreated: of the Q2 "recovered in PBMC" pathways,
#       which also move back toward HC in duodenum ("recovered in duodenum")
#       vs not, again direction-only relative to the ORIGINAL Q1 disease
#       direction.
#
# Each question is deliberately conditioned on the previous question's result
# set (a strict cascade / funnel), unlike the strict/directional sections in
# the Step 3 report, which each independently re-filter from the full pathway
# universe.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
})

# =============================================================================
# USER SETTINGS (mirrors qmd_step3_pathway_strict_then_directional.qmd)
# =============================================================================

pbmc_out_dir <- "deg_gsea_pooled_vs_within_pool"

duodenum_base_dir <- "/vf/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/OUTPUT_20251128_141815_WITH_BATCH/Monaco_fine"
duodenum_gsea_dir <- file.path(duodenum_base_dir, "gsea_results")

report_out <- file.path(pbmc_out_dir, "cascade_Mild_and_AllAPECED_vs_HC_recovery")
dir.create(report_out, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05  # only used for Q1 (defines "disease phenotype")

# Contrast names / folders
# Q1a now uses the direct subject-level pseudobulk outputs generated together
# with the All-APECED analysis.
pbmc_disease_folder    <- "all_apeced_vs_hc"
pbmc_disease_contrast  <- "MildUntreated_vs_HC_subject_level"  # Q1a: Mild vs HC

# Additional Q1 disease definition requested:
# subject-level All APECED (Mild untreated + Severe) vs HC
pbmc_all_disease_folder   <- "all_apeced_vs_hc"
pbmc_all_disease_contrast <- "AllAPECED_vs_HC"             # Q1b: All APECED vs HC

pbmc_treatment_folder   <- "paired_cross_pool"
pbmc_treatment_contrast <- "MildTreated_vs_MildUntreated_paired"  # Q2: pre vs post PBMC

duodenum_treatment_contrast <- "Treated_vs_Untreated"       # Q3: pre vs post duodenum

databases <- c("GO_BP", "Reactome")

# =============================================================================
# HELPERS (same conventions as the Step 3 QMD)
# =============================================================================

safe_name <- function(x) {
  x <- gsub("[/\\\\]", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(suppressWarnings(data.table::fread(path)), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  if ("reason" %in% names(x)) return(NULL)
  as.data.frame(x)
}

# Same as read_csv_safe(), but instead of silently returning NULL for a
# pipeline-generated placeholder file (e.g. "Fewer than 2 complete paired
# subjects"), it returns the reason string so it can be reported explicitly
# rather than lumped in with genuinely missing/broken files.
read_csv_with_skip_reason <- function(path) {
  if (!file.exists(path)) return(list(data = NULL, skip_reason = "File does not exist"))
  x <- tryCatch(suppressWarnings(data.table::fread(path)), error = function(e) NULL)
  if (is.null(x)) return(list(data = NULL, skip_reason = "fread() failed to parse file"))
  if (nrow(x) == 0) return(list(data = NULL, skip_reason = "File has zero rows"))
  if ("reason" %in% names(x)) {
    reason_text <- if (nrow(x) >= 1 && ncol(x) >= 1) as.character(x[[1]][1]) else "Unknown (reason column present but empty)"
    return(list(data = NULL, skip_reason = reason_text))
  }
  list(data = as.data.frame(x), skip_reason = NA_character_)
}

get_pbmc_celltypes <- function(folder = "pooled") {
  d <- file.path(pbmc_out_dir, folder, "gsea_by_celltype")
  if (!dir.exists(d)) stop("Missing folder: ", d)
  basename(list.dirs(d, recursive = FALSE, full.names = TRUE))
}

find_pbmc_gsea_file <- function(folder, celltype_safe, contrast, database) {
  file.path(
    pbmc_out_dir, folder, "gsea_by_celltype", celltype_safe,
    paste0(celltype_safe, "_", contrast, "_", database, "_results.csv")
  )
}

# Explicit overrides for Monaco fine labels where the duodenum export's
# cluster naming does not match a plain underscore->space conversion of the
# PBMC-side safe name. Add to this table as new mismatches are found (spot
# check with list_duodenum_celltypes.sh / diagnose_duodenum_parsing.R).
# Key = PBMC celltype label (post celltype_from_safe), Value = exact
# duodenum cluster label as it appears in the duodenum filenames.
duodenum_celltype_label_overrides <- c(
  "Th1 Th17 cells" = "Th1_Th17 cells"
)

find_duodenum_gsea_file <- function(celltype_label, database) {
  db_tag <- dplyr::case_when(
    database == "GO_BP" ~ "GO",
    database == "Reactome" ~ "REACTOME",
    TRUE ~ database
  )
  
  duodenum_label <- dplyr::if_else(
    celltype_label %in% names(duodenum_celltype_label_overrides),
    duodenum_celltype_label_overrides[celltype_label],
    celltype_label
  )
  
  file.path(
    duodenum_gsea_dir,
    paste0("Monaco_fine__cluster_", duodenum_label, "_GSEA_", db_tag, "_", duodenum_treatment_contrast, ".csv")
  )
}

celltype_from_safe <- function(x) stringr::str_replace_all(x, "_", " ")
num_col <- function(x) suppressWarnings(as.numeric(x))

sign_label <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "Missing",
    x > 0 ~ "Up",
    x < 0 ~ "Down",
    TRUE ~ "Zero"
  )
}

sign_opposite <- function(x, y) {
  !is.na(x) & !is.na(y) & sign(x) != sign(y) & sign(x) != 0 & sign(y) != 0
}

sign_same <- function(x, y) {
  !is.na(x) & !is.na(y) & sign(x) == sign(y) & sign(x) != 0 & sign(y) != 0
}

celltype_safe_vec <- get_pbmc_celltypes(pbmc_disease_folder)
celltype_key <- tibble(
  celltype_safe = celltype_safe_vec,
  celltype = celltype_from_safe(celltype_safe_vec)
)


bar_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

save_table <- function(df, fname) {
  data.table::fwrite(df, file.path(report_out, fname))
}

save_plot <- function(p, fname_stub, width = 12, height = 10) {
  ggsave(file.path(report_out, paste0(fname_stub, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(report_out, paste0(fname_stub, ".pdf")), p, width = width, height = height)
}

# =============================================================================
# Q0. File availability + skip-reason audit
#     Runs BEFORE Q1-Q3 so it's clear up front which celltype x database x
#     contrast combinations are genuinely untestable (and why, per the
#     pipeline's own "reason" placeholder files) versus a script-side bug.
#     Distinguishes: file missing, file empty, fread parse failure, and
#     pipeline-reported skip reasons (e.g. "Fewer than 2 complete paired
#     subjects") from each other.
# =============================================================================

audit_one_file <- function(celltype, celltype_safe, database, stage, file_path) {
  result <- read_csv_with_skip_reason(file_path)
  tibble(
    celltype = celltype,
    celltype_safe = celltype_safe,
    database = database,
    stage = stage,
    file_path = file_path,
    file_exists = file.exists(file_path),
    n_rows = if (is.null(result$data)) NA_integer_ else nrow(result$data),
    skip_reason = result$skip_reason
  )
}

file_availability_audit <- purrr::pmap_dfr(
  tidyr::expand_grid(celltype_key, database = databases),
  function(celltype, celltype_safe, database) {
    dplyr::bind_rows(
      audit_one_file(
        celltype, celltype_safe, database, "Q1a PBMC Mild disease",
        find_pbmc_gsea_file(
          pbmc_disease_folder,
          celltype_safe,
          pbmc_disease_contrast,
          database
        )
      ),
      audit_one_file(
        celltype, celltype_safe, database, "Q1b PBMC All APECED disease",
        find_pbmc_gsea_file(
          pbmc_all_disease_folder,
          celltype_safe,
          pbmc_all_disease_contrast,
          database
        )
      ),
      audit_one_file(
        celltype, celltype_safe, database, "Q2 PBMC treatment",
        find_pbmc_gsea_file(pbmc_treatment_folder, celltype_safe, pbmc_treatment_contrast, database)
      ),
      audit_one_file(
        celltype, celltype_safe, database, "Q3 Duodenum treatment",
        find_duodenum_gsea_file(celltype, database)
      )
    )
  }
)

save_table(file_availability_audit, "Q0_file_availability_and_skip_reasons.csv")

untestable_summary <- file_availability_audit %>%
  filter(!is.na(skip_reason)) %>%
  count(stage, skip_reason, celltype, database, name = "n_files") %>%
  arrange(stage, skip_reason, celltype, database)

save_table(untestable_summary, "Q0_untestable_celltypes_summary.csv")

cat("=== Q0 file availability audit ===\n")
cat("Total celltype x database x stage combinations checked:", nrow(file_availability_audit), "\n")
cat("Combinations with a skip reason (untestable):", sum(!is.na(file_availability_audit$skip_reason)), "\n\n")
if (nrow(untestable_summary) > 0) {
  cat("Untestable celltypes by stage and reason:\n")
  print(as.data.frame(untestable_summary))
}
cat("\n")


# =============================================================================
# Q1a. PBMC Mild untreated vs HC -- subject-level disease phenotype
# Q1b. PBMC All APECED vs HC -- requested expanded disease phenotype
#
# Both use FDR-significant pathways, split Up vs Down, per Monaco fine cell type.
# Positive NES means higher in the first group named in the contrast:
#   MildUntreated_vs_HC  -> higher in Mild untreated
#   AllAPECED_vs_HC      -> higher in All APECED
# =============================================================================

read_disease_phenotype_generic <- function(
    celltype_safe,
    celltype,
    database,
    folder,
    contrast,
    comparison_label
) {
  f <- find_pbmc_gsea_file(folder, celltype_safe, contrast, database)
  x <- read_csv_safe(f)
  if (is.null(x)) return(NULL)
  
  needed <- c("ID", "Description", "NES", "p.adjust")
  if (!all(needed %in% names(x))) return(NULL)
  
  x %>%
    transmute(
      celltype = celltype,
      celltype_safe = celltype_safe,
      database = database,
      comparison = comparison_label,
      feature_id = ID,
      feature_label = Description,
      disease_NES = num_col(NES),
      disease_FDR = num_col(p.adjust),
      disease_file = f
    ) %>%
    filter(
      !is.na(disease_FDR),
      disease_FDR < fdr_cutoff,
      !is.na(disease_NES),
      disease_NES != 0
    ) %>%
    mutate(disease_direction = sign_label(disease_NES))
}

read_mild_disease_phenotype <- function(celltype_safe, celltype, database) {
  read_disease_phenotype_generic(
    celltype_safe = celltype_safe,
    celltype = celltype,
    database = database,
    folder = pbmc_disease_folder,
    contrast = pbmc_disease_contrast,
    comparison_label = "Mild untreated vs HC"
  )
}

read_all_apeced_disease_phenotype <- function(celltype_safe, celltype, database) {
  read_disease_phenotype_generic(
    celltype_safe = celltype_safe,
    celltype = celltype,
    database = database,
    folder = pbmc_all_disease_folder,
    contrast = pbmc_all_disease_contrast,
    comparison_label = "All APECED vs HC"
  )
}

# ---------------------------------------------------------------------------
# Q1a: Mild untreated vs HC
# ---------------------------------------------------------------------------
disease_phenotype <- purrr::pmap_dfr(
  tidyr::expand_grid(celltype_key, database = databases),
  read_mild_disease_phenotype
)

if (nrow(disease_phenotype) == 0) {
  warning(
    "Q1a: no Mild-vs-HC disease-phenotype pathways were parsed. ",
    "Check pbmc_disease_contrast and file paths."
  )
}

save_table(
  disease_phenotype,
  "Q1a_disease_phenotype_pathways_MildUntreated_vs_HC.csv"
)

q1_counts <- disease_phenotype %>%
  count(database, celltype, disease_direction, name = "n") %>%
  mutate(
    disease_direction = factor(
      disease_direction,
      levels = c("Up", "Down")
    )
  ) %>%
  complete(
    database,
    celltype,
    disease_direction,
    fill = list(n = 0)
  )

save_table(
  q1_counts,
  "Q1a_disease_phenotype_counts_MildUntreated_vs_HC.csv"
)

q1_celltype_order <- q1_counts %>%
  group_by(celltype) %>%
  summarise(total = sum(n), .groups = "drop")

q1_plot_df <- q1_counts %>%
  left_join(q1_celltype_order, by = "celltype") %>%
  mutate(celltype = forcats::fct_reorder(celltype, total))

p1 <- ggplot(
  q1_plot_df,
  aes(x = n, y = celltype, fill = disease_direction)
) +
  geom_col(position = "stack", width = 0.75) +
  facet_wrap(~ database, nrow = 1) +
  scale_fill_manual(values = c(Up = "#d73027", Down = "#4575b4")) +
  labs(
    title = "Q1a: PBMC disease phenotype — Mild untreated vs HC",
    subtitle = "FDR < 0.05; significant pathways per Monaco fine population",
    x = "Number of significant pathways",
    y = NULL,
    fill = "Direction in Mild vs HC"
  ) +
  bar_theme

save_plot(
  p1,
  "Q1a_disease_phenotype_MildUntreated_vs_HC_barplot"
)

# ---------------------------------------------------------------------------
# Q1b: All APECED vs HC
# ---------------------------------------------------------------------------
all_apeced_disease_phenotype <- purrr::pmap_dfr(
  tidyr::expand_grid(celltype_key, database = databases),
  read_all_apeced_disease_phenotype
)

if (nrow(all_apeced_disease_phenotype) == 0) {
  warning(
    "Q1b: no All-APECED-vs-HC disease-phenotype pathways were parsed. ",
    "Expected files under: ",
    file.path(
      pbmc_out_dir,
      pbmc_all_disease_folder,
      "gsea_by_celltype"
    )
  )
}

save_table(
  all_apeced_disease_phenotype,
  "Q1b_disease_phenotype_pathways_AllAPECED_vs_HC.csv"
)

q1_all_counts <- all_apeced_disease_phenotype %>%
  count(database, celltype, disease_direction, name = "n") %>%
  mutate(
    disease_direction = factor(
      disease_direction,
      levels = c("Up", "Down")
    )
  ) %>%
  complete(
    database,
    celltype,
    disease_direction,
    fill = list(n = 0)
  )

save_table(
  q1_all_counts,
  "Q1b_disease_phenotype_counts_AllAPECED_vs_HC.csv"
)

q1_all_celltype_order <- q1_all_counts %>%
  group_by(celltype) %>%
  summarise(total = sum(n), .groups = "drop")

q1_all_plot_df <- q1_all_counts %>%
  left_join(q1_all_celltype_order, by = "celltype") %>%
  mutate(celltype = forcats::fct_reorder(celltype, total))

p1_all <- ggplot(
  q1_all_plot_df,
  aes(x = n, y = celltype, fill = disease_direction)
) +
  geom_col(position = "stack", width = 0.75) +
  facet_wrap(~ database, nrow = 1) +
  scale_fill_manual(values = c(Up = "#d73027", Down = "#4575b4")) +
  labs(
    title = "Q1b: PBMC disease phenotype — All APECED vs HC",
    subtitle = "Subject-level pseudobulk; FDR < 0.05",
    x = "Number of significant pathways",
    y = NULL,
    fill = "Direction in All APECED vs HC"
  ) +
  bar_theme

save_plot(
  p1_all,
  "Q1b_disease_phenotype_AllAPECED_vs_HC_barplot"
)

# ---------------------------------------------------------------------------
# Q1 side-by-side comparison
# ---------------------------------------------------------------------------
q1_comparison_counts <- bind_rows(
  q1_counts %>%
    mutate(comparison = "Mild untreated vs HC"),
  q1_all_counts %>%
    mutate(comparison = "All APECED vs HC")
) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "Mild untreated vs HC",
        "All APECED vs HC"
      )
    )
  )

save_table(
  q1_comparison_counts,
  "Q1_side_by_side_disease_phenotype_counts.csv"
)

q1_comparison_order <- q1_comparison_counts %>%
  group_by(celltype) %>%
  summarise(total = sum(n), .groups = "drop")

q1_comparison_plot_df <- q1_comparison_counts %>%
  left_join(q1_comparison_order, by = "celltype") %>%
  mutate(celltype = forcats::fct_reorder(celltype, total))

p1_side_by_side <- ggplot(
  q1_comparison_plot_df,
  aes(x = n, y = celltype, fill = disease_direction)
) +
  geom_col(position = "stack", width = 0.75) +
  facet_grid(
    comparison ~ database,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(values = c(Up = "#d73027", Down = "#4575b4")) +
  labs(
    title = "Q1: PBMC disease-phenotype pathway counts",
    subtitle = "Mild untreated vs HC compared with All APECED vs HC; FDR < 0.05",
    x = "Number of significant pathways",
    y = NULL,
    fill = "Disease direction"
  ) +
  bar_theme

save_plot(
  p1_side_by_side,
  "Q1_side_by_side_Mild_vs_AllAPECED_disease_phenotype",
  width = 14,
  height = 16
)

# Pathway-level concordance table between the two Q1 definitions
q1_pathway_comparison <- full_join(
  disease_phenotype %>%
    select(
      celltype,
      celltype_safe,
      database,
      feature_id,
      feature_label,
      mild_NES = disease_NES,
      mild_FDR = disease_FDR,
      mild_direction = disease_direction
    ),
  all_apeced_disease_phenotype %>%
    select(
      celltype,
      celltype_safe,
      database,
      feature_id,
      feature_label,
      all_apeced_NES = disease_NES,
      all_apeced_FDR = disease_FDR,
      all_apeced_direction = disease_direction
    ),
  by = c(
    "celltype",
    "celltype_safe",
    "database",
    "feature_id",
    "feature_label"
  )
) %>%
  mutate(
    significant_in_mild = !is.na(mild_FDR),
    significant_in_all_apeced = !is.na(all_apeced_FDR),
    significance_overlap = case_when(
      significant_in_mild & significant_in_all_apeced ~ "Significant in both",
      significant_in_mild & !significant_in_all_apeced ~ "Mild only",
      !significant_in_mild & significant_in_all_apeced ~ "All APECED only",
      TRUE ~ "Neither"
    ),
    direction_concordant = case_when(
      !is.na(mild_NES) & !is.na(all_apeced_NES) ~ sign(mild_NES) == sign(all_apeced_NES),
      TRUE ~ NA
    )
  )

save_table(
  q1_pathway_comparison,
  "Q1_pathway_level_Mild_vs_AllAPECED_comparison.csv"
)

# =============================================================================
# Q2. PBMC MildTreated vs MildUntreated -- "recovered in PBMC"
#     NOTE: Q2/Q3 continue to use Q1a Mild-vs-HC as the original disease
#     phenotype, preserving the original cascade logic. Q1b is an added
#     robustness comparison and does not replace Q1a in the recovery cascade.
#     Restricted to the Q1 disease-phenotype pathways (same celltype+database
#     +feature_id). Direction only, no significance filter on the treatment
#     stat. "Recovered" = treatment NES sign opposite to the ORIGINAL disease
#     NES sign (i.e. moving back toward HC).
# =============================================================================

read_pbmc_treatment_join <- function(celltype_safe, celltype, database) {
  f <- find_pbmc_gsea_file(pbmc_treatment_folder, celltype_safe, pbmc_treatment_contrast, database)
  x <- read_csv_safe(f)
  if (is.null(x)) return(NULL)
  
  needed <- c("ID", "NES")
  if (!all(needed %in% names(x))) return(NULL)
  
  x %>%
    transmute(
      celltype = celltype,
      celltype_safe = celltype_safe,
      database = database,
      feature_id = ID,
      pbmc_treatment_NES = num_col(NES),
      pbmc_treatment_file = f
    )
}

pbmc_treatment_all <- purrr::pmap_dfr(
  tidyr::expand_grid(celltype_key, database = databases),
  read_pbmc_treatment_join
)

pbmc_recovery <- disease_phenotype %>%
  left_join(
    pbmc_treatment_all,
    by = c("celltype", "celltype_safe", "database", "feature_id")
  ) %>%
  mutate(
    tested_in_pbmc_treatment = !is.na(pbmc_treatment_NES),
    pbmc_recovery_status = dplyr::case_when(
      !tested_in_pbmc_treatment ~ "Not tested in PBMC treatment",
      sign_opposite(disease_NES, pbmc_treatment_NES) ~ "Recovered after treatment in PBMC",
      sign_same(disease_NES, pbmc_treatment_NES) ~ "Not recovered after treatment in PBMC",
      TRUE ~ "Other"
    )
  )

pbmc_treatment_skip_reasons <- file_availability_audit %>%
  filter(stage == "Q2 PBMC treatment") %>%
  select(celltype, database, pbmc_treatment_skip_reason = skip_reason)

pbmc_recovery <- pbmc_recovery %>%
  left_join(pbmc_treatment_skip_reasons, by = c("celltype", "database"))

save_table(pbmc_recovery, "Q2_pbmc_treatment_recovery_status.csv")

q2_counts <- pbmc_recovery %>%
  count(database, celltype, pbmc_recovery_status, name = "n") %>%
  mutate(
    pbmc_recovery_status = factor(
      pbmc_recovery_status,
      levels = c(
        "Recovered after treatment in PBMC",
        "Not recovered after treatment in PBMC",
        "Not tested in PBMC treatment"
      )
    )
  ) %>%
  complete(database, celltype, pbmc_recovery_status, fill = list(n = 0))

save_table(q2_counts, "Q2_pbmc_treatment_recovery_counts_per_celltype.csv")

q2_celltype_order <- q2_counts %>%
  group_by(celltype) %>%
  summarise(total = sum(n), .groups = "drop")

q2_plot_df <- q2_counts %>%
  left_join(q2_celltype_order, by = "celltype") %>%
  mutate(celltype = forcats::fct_reorder(celltype, total))

p2 <- ggplot(q2_plot_df, aes(x = n, y = celltype, fill = pbmc_recovery_status)) +
  geom_col(position = "stack", width = 0.75) +
  facet_wrap(~ database, nrow = 1) +
  scale_fill_manual(values = c(
    "Recovered after treatment in PBMC" = "#1a9850",
    "Not recovered after treatment in PBMC" = "#d73027",
    "Not tested in PBMC treatment" = "grey70"
  )) +
  labs(
    title = "Q2: PBMC pre vs post treatment recovery of disease-phenotype pathways",
    subtitle = "Restricted to Q1 disease-phenotype pathways; direction only (no significance filter)",
    x = "Number of pathways", y = NULL, fill = NULL
  ) +
  bar_theme

save_plot(p2, "Q2_pbmc_treatment_recovery_barplot")

# =============================================================================
# Q3. Duodenum Treated vs Untreated -- "recovered in duodenum"
#     Restricted to the Q2 "Recovered after treatment in PBMC" pathways.
#     Direction only, relative to the ORIGINAL Q1 disease NES sign.
# =============================================================================

read_duodenum_treatment_join <- function(celltype, database) {
  f <- find_duodenum_gsea_file(celltype, database)
  x <- read_csv_safe(f)
  if (is.null(x)) return(NULL)
  
  needed <- c("ID", "NES")
  if (!all(needed %in% names(x))) return(NULL)
  
  x %>%
    transmute(
      celltype = celltype,
      database = database,
      feature_id = ID,
      duodenum_treatment_NES = num_col(NES),
      duodenum_treatment_file = f
    )
}

pbmc_recovered_only <- pbmc_recovery %>%
  filter(pbmc_recovery_status == "Recovered after treatment in PBMC")

duodenum_treatment_all <- purrr::pmap_dfr(
  tidyr::expand_grid(distinct(celltype_key, celltype), database = databases),
  read_duodenum_treatment_join
)

duodenum_recovery <- pbmc_recovered_only %>%
  left_join(
    duodenum_treatment_all,
    by = c("celltype", "database", "feature_id")
  ) %>%
  mutate(
    tested_in_duodenum_treatment = !is.na(duodenum_treatment_NES),
    duodenum_recovery_status = dplyr::case_when(
      !tested_in_duodenum_treatment ~ "Not tested in duodenum treatment",
      sign_opposite(disease_NES, duodenum_treatment_NES) ~ "Recovered after treatment in Duodenum",
      sign_same(disease_NES, duodenum_treatment_NES) ~ "Not recovered after treatment in Duodenum",
      TRUE ~ "Other"
    )
  )

duodenum_treatment_skip_reasons <- file_availability_audit %>%
  filter(stage == "Q3 Duodenum treatment") %>%
  select(celltype, database, duodenum_treatment_skip_reason = skip_reason)

duodenum_recovery <- duodenum_recovery %>%
  left_join(duodenum_treatment_skip_reasons, by = c("celltype", "database"))

save_table(duodenum_recovery, "Q3_duodenum_treatment_recovery_status.csv")

q3_counts <- duodenum_recovery %>%
  count(database, celltype, duodenum_recovery_status, name = "n") %>%
  mutate(
    duodenum_recovery_status = factor(
      duodenum_recovery_status,
      levels = c(
        "Recovered after treatment in Duodenum",
        "Not recovered after treatment in Duodenum",
        "Not tested in duodenum treatment"
      )
    )
  ) %>%
  complete(database, celltype, duodenum_recovery_status, fill = list(n = 0))

save_table(q3_counts, "Q3_duodenum_treatment_recovery_counts_per_celltype.csv")

q3_celltype_order <- q3_counts %>%
  group_by(celltype) %>%
  summarise(total = sum(n), .groups = "drop")

q3_plot_df <- q3_counts %>%
  left_join(q3_celltype_order, by = "celltype") %>%
  mutate(celltype = forcats::fct_reorder(celltype, total))

p3 <- ggplot(q3_plot_df, aes(x = n, y = celltype, fill = duodenum_recovery_status)) +
  geom_col(position = "stack", width = 0.75) +
  facet_wrap(~ database, nrow = 1) +
  scale_fill_manual(values = c(
    "Recovered after treatment in Duodenum" = "#1a9850",
    "Not recovered after treatment in Duodenum" = "#d73027",
    "Not tested in duodenum treatment" = "grey70"
  )) +
  labs(
    title = "Q3: Duodenum pre vs post treatment recovery",
    subtitle = "Restricted to Q2 'Recovered after treatment in PBMC' pathways; direction relative to original Q1 disease effect",
    x = "Number of pathways", y = NULL, fill = NULL
  ) +
  bar_theme

save_plot(p3, "Q3_duodenum_treatment_recovery_barplot")

# =============================================================================
# Console summary
# =============================================================================

cat("\n=== Cascade recovery analysis complete ===\n")
cat("Output directory:", normalizePath(report_out, mustWork = FALSE), "\n\n")

cat(
  "Q1a disease-phenotype pathways (Mild untreated vs HC, FDR <",
  fdr_cutoff, "):", nrow(disease_phenotype), "\n"
)
cat(
  "Q1b disease-phenotype pathways (All APECED vs HC, FDR <",
  fdr_cutoff, "):", nrow(all_apeced_disease_phenotype), "\n"
)
cat("Q2 pathways tested in PBMC treatment:", sum(pbmc_recovery$tested_in_pbmc_treatment), "\n")
cat("Q2 'Recovered after treatment in PBMC':",
    sum(pbmc_recovery$pbmc_recovery_status == "Recovered after treatment in PBMC"), "\n")
cat("Q3 pathways tested in duodenum treatment:", sum(duodenum_recovery$tested_in_duodenum_treatment), "\n")
cat("Q3 'Recovered after treatment in Duodenum':",
    sum(duodenum_recovery$duodenum_recovery_status == "Recovered after treatment in Duodenum", na.rm = TRUE), "\n")

cat("\nFiles written:\n")
print(list.files(report_out))