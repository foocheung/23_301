# =============================================================================
# YU2021 HEALTHY + INTERNAL UNTREATED/TREATED
# COMBINED PSEUDOBULK LIMMA-VOOM + FGSEA PIPELINE
#
# CONTRASTS
#   1. Untreated_vs_Healthy
#        Unpaired cross-study comparison
#
#   2. Treated_vs_Healthy
#        Unpaired cross-study comparison
#
#   3. Treated_vs_Untreated
#        Explicitly paired internal comparison:
#        design = ~ Subject_fixed + Treatment
#
# BIOLOGICAL REPLICATES
#   Yu2021   = donorID_unified
#   Internal = Subject_fixed
#
# CELL-TYPE ANNOTATIONS
#   Yu2021   = level_2_annot
#   Internal = PanGI_L2
#
# PATHWAYS
#   Hallmark
#   Reactome
#   GO Biological Process
#
# IMPORTANT
#   - Uses original RNA counts.
#   - Does not use Harmony values for DEG.
#   - Healthy comparisons remain confounded with study.
#   - Treated_vs_Untreated is a true paired subject-level contrast.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(scales)
  library(fs)
})

set.seed(1234)

# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG <- list(
  
  # ---------------------------------------------------------------------------
  # Input
  # ---------------------------------------------------------------------------
  
  yu_rds = file.path(
    "YU2021_INTERNAL_HARMONY",
    "Yu2021_independent_processed.rds"
  ),
  
  internal_rds = file.path(
    "YU2021_INTERNAL_HARMONY",
    "Internal_independent_processed.rds"
  ),
  
  assay = "RNA",
  
  # ---------------------------------------------------------------------------
  # Output
  # ---------------------------------------------------------------------------
  
  out_dir =
    "deg_gsea_YuHealthy_Internal_paired_limma_fgsea",
  
  # ---------------------------------------------------------------------------
  # Pseudobulk requirements
  # ---------------------------------------------------------------------------
  
  min_cells_per_pseudobulk = 20,
  
  min_subjects_per_group = 2,
  
  min_paired_subjects = 2,
  
  # ---------------------------------------------------------------------------
  # Gene filtering
  # ---------------------------------------------------------------------------
  
  min_count = 10,
  
  min_total_count = 20,
  
  min_genes_after_filtering = 100,
  
  # ---------------------------------------------------------------------------
  # Limma
  # ---------------------------------------------------------------------------
  
  robust_ebayes = TRUE,
  
  deg_fdr = 0.10,
  
  top_n_gene_labels = 10,
  
  # ---------------------------------------------------------------------------
  # FGSEA
  # ---------------------------------------------------------------------------
  
  gsea_min_size = 15,
  
  gsea_max_size = 500,
  
  gsea_eps = 0,
  
  pathway_fdr = 0.10,
  
  strict_pathway_fdr = 0.05,
  
  top_n_pathways = 10,
  
  # ---------------------------------------------------------------------------
  # Output options
  # ---------------------------------------------------------------------------
  
  save_pseudobulk_counts = TRUE,
  
  save_model_objects = FALSE,
  
  make_volcano_plots = TRUE,
  
  make_pathway_plots = TRUE
)

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================

DEG_DIR <- file.path(
  CONFIG$out_dir,
  "01_DEG"
)

PB_DIR <- file.path(
  CONFIG$out_dir,
  "02_PSEUDOBULK"
)

GSEA_DIR <- file.path(
  CONFIG$out_dir,
  "03_GSEA"
)

SUMMARY_DIR <- file.path(
  CONFIG$out_dir,
  "04_SUMMARY"
)

FIGURE_DIR <- file.path(
  CONFIG$out_dir,
  "05_FIGURES"
)

DECISION_DIR <- file.path(
  SUMMARY_DIR,
  "PATHWAY_DECISION_TREE"
)

dir_create(DEG_DIR)
dir_create(PB_DIR)
dir_create(GSEA_DIR)
dir_create(SUMMARY_DIR)
dir_create(FIGURE_DIR)
dir_create(DECISION_DIR)

# =============================================================================
# LOGGING
# =============================================================================

logmsg <- function(fmt, ...) {
  
  message(
    sprintf(
      paste0(
        "[",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "] ",
        fmt
      ),
      ...
    )
  )
}

# =============================================================================
# GENERAL HELPERS
# =============================================================================

clean_character <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x[
    is.na(x) |
      x == "" |
      tolower(x) %in% c(
        "na",
        "n/a",
        "none",
        "null",
        "missing",
        "unknown"
      )
  ] <- NA_character_
  
  x
}

clean_name <- function(x) {
  
  x %>%
    as.character() %>%
    str_replace_all(
      "[^A-Za-z0-9_.-]",
      "_"
    ) %>%
    str_replace_all(
      "_+",
      "_"
    ) %>%
    str_replace(
      "^_",
      ""
    ) %>%
    str_replace(
      "_$",
      ""
    )
}

save_plot_both <- function(
    plot_object,
    filename_prefix,
    width,
    height
) {
  
  ggsave(
    filename = paste0(
      filename_prefix,
      ".png"
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(
      filename_prefix,
      ".pdf"
    ),
    plot = plot_object,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}

# =============================================================================
# RETRIEVE COUNTS SAFELY
# =============================================================================

get_counts_matrix <- function(object) {
  
  DefaultAssay(object) <- CONFIG$assay
  
  if (
    inherits(
      object[[CONFIG$assay]],
      "Assay5"
    )
  ) {
    
    count_layers <- grep(
      "^counts",
      Layers(
        object[[CONFIG$assay]]
      ),
      value = TRUE
    )
    
    if (length(count_layers) > 1) {
      
      object <- JoinLayers(
        object,
        assay = CONFIG$assay
      )
    }
  }
  
  counts <- tryCatch(
    {
      GetAssayData(
        object,
        assay = CONFIG$assay,
        layer = "counts"
      )
    },
    error = function(e) {
      
      GetAssayData(
        object,
        assay = CONFIG$assay,
        slot = "counts"
      )
    }
  )
  
  list(
    object = object,
    counts = counts
  )
}

# =============================================================================
# COLLAPSE CELLS INTO PSEUDOBULKS
# =============================================================================

collapse_sparse_columns <- function(
    counts,
    groups
) {
  
  groups <- factor(groups)
  
  design <- sparse.model.matrix(
    ~ 0 + groups
  )
  
  colnames(design) <- levels(groups)
  
  counts %*% design
}

# =============================================================================
# LOAD OBJECTS
# =============================================================================

logmsg(
  "Loading Yu2021 object"
)

if (!file.exists(CONFIG$yu_rds)) {
  stop(
    "Yu RDS not found: ",
    CONFIG$yu_rds
  )
}

#yu <- readRDS(
#  CONFIG$yu_rds
#)

logmsg(
  "Loading internal object"
)

if (!file.exists(CONFIG$internal_rds)) {
  stop(
    "Internal RDS not found: ",
    CONFIG$internal_rds
  )
}

#obj <- readRDS(
#  CONFIG$internal_rds
#)

if (!inherits(yu, "Seurat")) {
  stop("Yu input is not a Seurat object.")
}

if (!inherits(obj, "Seurat")) {
  stop("Internal input is not a Seurat object.")
}

DefaultAssay(yu) <- CONFIG$assay
DefaultAssay(obj) <- CONFIG$assay

# =============================================================================
# VERIFY REQUIRED METADATA
# =============================================================================

required_yu <- c(
  "sampleID",
  "donorID_unified",
  "level_2_annot"
)

required_internal <- c(
  "orig.ident",
  "Subject_fixed",
  "Treatment",
  "PanGI_L2"
)

missing_yu <- setdiff(
  required_yu,
  colnames(yu@meta.data)
)

missing_internal <- setdiff(
  required_internal,
  colnames(obj@meta.data)
)

if (length(missing_yu) > 0) {
  stop(
    "Yu object is missing: ",
    paste(
      missing_yu,
      collapse = ", "
    )
  )
}

if (length(missing_internal) > 0) {
  stop(
    "Internal object is missing: ",
    paste(
      missing_internal,
      collapse = ", "
    )
  )
}

# =============================================================================
# STANDARDIZE METADATA
# =============================================================================

yu$dataset <- "Yu2021"

yu$technical_sample <- clean_character(
  yu$sampleID
)

yu$biological_subject <- clean_character(
  yu$donorID_unified
)

yu$condition_common <- "Healthy"

yu$celltype_common <- clean_character(
  yu$level_2_annot
)

yu$subject_uid <- paste(
  "Yu2021",
  yu$biological_subject,
  sep = "__"
)

obj$dataset <- "Internal"

obj$technical_sample <- clean_character(
  obj$orig.ident
)

obj$biological_subject <- clean_character(
  obj$Subject_fixed
)

obj$condition_common <- clean_character(
  obj$Treatment
)

obj$celltype_common <- clean_character(
  obj$PanGI_L2
)

obj$subject_uid <- paste(
  "Internal",
  obj$biological_subject,
  sep = "__"
)

# Keep internal Treated and Untreated only
obj <- subset(
  obj,
  subset =
    condition_common %in%
    c(
      "Untreated",
      "Treated"
    )
)

# =============================================================================
# MAKE CELL NAMES UNIQUE
# =============================================================================

if (!all(
  startsWith(
    colnames(yu),
    "YU2021_"
  )
)) {
  
  yu <- RenameCells(
    yu,
    add.cell.id = "YU2021"
  )
}

if (!all(
  startsWith(
    colnames(obj),
    "INTERNAL_"
  )
)) {
  
  obj <- RenameCells(
    obj,
    add.cell.id = "INTERNAL"
  )
}

if (
  length(
    intersect(
      colnames(yu),
      colnames(obj)
    )
  ) > 0
) {
  
  stop(
    "Cell barcodes overlap after prefixing."
  )
}

# =============================================================================
# RESTRICT TO COMMON GENES AND MERGE
# =============================================================================

common_genes <- intersect(
  rownames(yu),
  rownames(obj)
)

logmsg(
  "Shared genes: %d",
  length(common_genes)
)

if (length(common_genes) < 1000) {
  stop(
    "Too few shared genes: ",
    length(common_genes)
  )
}

yu <- subset(
  yu,
  features = common_genes
)

obj <- subset(
  obj,
  features = common_genes
)

combined <- merge(
  x = yu,
  y = obj,
  merge.data = FALSE,
  project = "Yu_Internal_Pseudobulk"
)

DefaultAssay(combined) <- CONFIG$assay

count_result <- get_counts_matrix(
  combined
)

combined <- count_result$object
counts_matrix <- count_result$counts

# =============================================================================
# ALIGN COUNTS AND METADATA
# =============================================================================

cell_order <- intersect(
  colnames(counts_matrix),
  rownames(combined@meta.data)
)

counts_matrix <- counts_matrix[
  ,
  cell_order,
  drop = FALSE
]

metadata <- combined@meta.data[
  cell_order,
  ,
  drop = FALSE
] %>%
  rownames_to_column(
    "cell_barcode"
  ) %>%
  transmute(
    cell_barcode,
    dataset =
      clean_character(dataset),
    technical_sample =
      clean_character(technical_sample),
    biological_subject =
      clean_character(biological_subject),
    subject_uid =
      clean_character(subject_uid),
    condition =
      clean_character(condition_common),
    celltype =
      clean_character(celltype_common)
  ) %>%
  filter(
    !is.na(dataset),
    !is.na(subject_uid),
    !is.na(condition),
    !is.na(celltype),
    condition %in%
      c(
        "Healthy",
        "Untreated",
        "Treated"
      )
  ) %>%
  mutate(
    celltype_clean =
      clean_name(celltype),
    
    pseudobulk_id = paste(
      celltype_clean,
      subject_uid,
      condition,
      sep = "__"
    )
  )

counts_matrix <- counts_matrix[
  ,
  metadata$cell_barcode,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts_matrix),
    metadata$cell_barcode
  )
)

# =============================================================================
# DESIGN OVERVIEW
# =============================================================================

design_overview <- metadata %>%
  distinct(
    dataset,
    subject_uid,
    biological_subject,
    condition
  ) %>%
  count(
    dataset,
    condition,
    name = "n_subjects"
  )

print(
  design_overview
)

write_csv(
  design_overview,
  file.path(
    SUMMARY_DIR,
    "subject_design_overview.csv"
  )
)

cat("\nInternal subject pairing:\n")

print(
  metadata %>%
    filter(
      dataset == "Internal"
    ) %>%
    distinct(
      biological_subject,
      condition
    ) %>%
    table()
)

# =============================================================================
# CELL COUNTS PER PSEUDOBULK
# =============================================================================

cell_count_summary <- metadata %>%
  count(
    celltype,
    celltype_clean,
    dataset,
    subject_uid,
    biological_subject,
    condition,
    pseudobulk_id,
    name = "n_cells"
  )

write_csv(
  cell_count_summary,
  file.path(
    SUMMARY_DIR,
    "cells_per_pseudobulk.csv"
  )
)

# =============================================================================
# VOLCANO PLOT
# =============================================================================

make_volcano <- function(
    result,
    contrast,
    celltype
) {
  
  plot_data <- result %>%
    mutate(
      significance = case_when(
        
        adj.P.Val <= CONFIG$deg_fdr &
          logFC > 0 ~
          "Up",
        
        adj.P.Val <= CONFIG$deg_fdr &
          logFC < 0 ~
          "Down",
        
        TRUE ~
          "Not significant"
      ),
      
      minus_log10_fdr =
        -log10(
          pmax(
            adj.P.Val,
            1e-300
          )
        )
    )
  
  label_data <- plot_data %>%
    filter(
      significance !=
        "Not significant"
    ) %>%
    arrange(
      adj.P.Val,
      desc(abs(logFC))
    ) %>%
    slice_head(
      n =
        CONFIG$top_n_gene_labels
    )
  
  ggplot(
    plot_data,
    aes(
      x = logFC,
      y = minus_log10_fdr
    )
  ) +
    geom_point(
      aes(
        color = significance
      ),
      size = 0.8,
      alpha = 0.7
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept =
        -log10(CONFIG$deg_fdr),
      linetype = "dashed"
    ) +
    geom_text_repel(
      data = label_data,
      aes(
        label = gene
      ),
      size = 3,
      max.overlaps = 30
    ) +
    scale_color_manual(
      values = c(
        Up = "#C0392B",
        Down = "#2980B9",
        `Not significant` =
          "grey75"
      )
    ) +
    labs(
      title = paste0(
        celltype,
        ": ",
        contrast
      ),
      subtitle =
        "Pseudobulk limma-voom",
      x =
        "log2 fold change",
      y =
        expression(
          -log[10](FDR)
        ),
      color = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      plot.title =
        element_text(
          face = "bold"
        ),
      panel.grid.minor =
        element_blank(),
      legend.position =
        "bottom"
    )
}

# =============================================================================
# FIT AN UNPAIRED TWO-GROUP MODEL
#
# Used for:
#   Untreated_vs_Healthy
#   Treated_vs_Healthy
# =============================================================================

fit_unpaired_contrast <- function(
    pb_counts,
    pb_metadata,
    numerator,
    denominator,
    celltype_value,
    celltype_clean
) {
  
  contrast_name <- paste0(
    numerator,
    "_vs_",
    denominator
  )
  
  current_metadata <- pb_metadata %>%
    filter(
      condition %in%
        c(
          denominator,
          numerator
        )
    ) %>%
    mutate(
      condition = factor(
        condition,
        levels = c(
          denominator,
          numerator
        )
      )
    )
  
  subject_counts <- current_metadata %>%
    count(
      condition,
      name = "n_subjects"
    )
  
  n_numerator <- sum(
    current_metadata$condition ==
      numerator
  )
  
  n_denominator <- sum(
    current_metadata$condition ==
      denominator
  )
  
  if (
    n_numerator <
    CONFIG$min_subjects_per_group ||
    n_denominator <
    CONFIG$min_subjects_per_group
  ) {
    
    return(
      list(
        status = "Skipped",
        reason = paste0(
          "Insufficient subjects: ",
          numerator,
          "=",
          n_numerator,
          "; ",
          denominator,
          "=",
          n_denominator
        ),
        result = NULL
      )
    )
  }
  
  current_counts <- pb_counts[
    ,
    current_metadata$pseudobulk_id,
    drop = FALSE
  ]
  
  design <- model.matrix(
    ~ 0 + condition,
    data = current_metadata
  )
  
  colnames(design) <- levels(
    current_metadata$condition
  )
  
  dge <- DGEList(
    counts = current_counts
  )
  
  keep <- filterByExpr(
    dge,
    design = design,
    min.count =
      CONFIG$min_count,
    min.total.count =
      CONFIG$min_total_count
  )
  
  if (
    sum(keep) <
    CONFIG$min_genes_after_filtering
  ) {
    
    return(
      list(
        status = "Skipped",
        reason = paste0(
          "Only ",
          sum(keep),
          " genes passed filtering"
        ),
        result = NULL
      )
    )
  }
  
  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  dge <- calcNormFactors(
    dge,
    method = "TMM"
  )
  
  voom_object <- voom(
    dge,
    design,
    plot = FALSE
  )
  
  fit <- lmFit(
    voom_object,
    design
  )
  
  contrast_matrix <- makeContrasts(
    contrasts = paste0(
      numerator,
      "-",
      denominator
    ),
    levels = design
  )
  
  fit <- contrasts.fit(
    fit,
    contrast_matrix
  )
  
  fit <- eBayes(
    fit,
    robust =
      CONFIG$robust_ebayes
  )
  
  result <- topTable(
    fit,
    coef = 1,
    number = Inf,
    sort.by = "P"
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    mutate(
      contrast = contrast_name,
      numerator = numerator,
      denominator = denominator,
      paired = FALSE,
      celltype = celltype_value,
      celltype_clean =
        celltype_clean,
      n_numerator_subjects =
        n_numerator,
      n_denominator_subjects =
        n_denominator
    )
  
  list(
    status = "Completed",
    reason = NA_character_,
    result = result,
    voom = voom_object,
    fit = fit,
    design = design
  )
}

# =============================================================================
# FIT EXPLICITLY PAIRED TREATED VS UNTREATED MODEL
#
# Internal subjects only.
#
# Model:
#   ~ biological_subject + condition
#
# condition coefficient:
#   Treated relative to Untreated
# =============================================================================

fit_paired_treatment_contrast <- function(
    pb_counts,
    pb_metadata,
    celltype_value,
    celltype_clean
) {
  
  paired_metadata <- pb_metadata %>%
    filter(
      dataset == "Internal",
      condition %in%
        c(
          "Untreated",
          "Treated"
        )
    ) %>%
    group_by(
      biological_subject
    ) %>%
    filter(
      n_distinct(condition) == 2
    ) %>%
    ungroup()
  
  paired_subjects <- sort(
    unique(
      paired_metadata$biological_subject
    )
  )
  
  n_pairs <- length(
    paired_subjects
  )
  
  if (
    n_pairs <
    CONFIG$min_paired_subjects
  ) {
    
    return(
      list(
        status = "Skipped",
        reason = paste0(
          "Only ",
          n_pairs,
          " complete paired subjects"
        ),
        result = NULL
      )
    )
  }
  
  paired_metadata <- paired_metadata %>%
    mutate(
      biological_subject =
        factor(
          biological_subject
        ),
      
      condition =
        factor(
          condition,
          levels = c(
            "Untreated",
            "Treated"
          )
        )
    ) %>%
    arrange(
      biological_subject,
      condition
    )
  
  # Confirm one untreated and one treated sample per retained subject
  pair_table <- table(
    paired_metadata$biological_subject,
    paired_metadata$condition
  )
  
  valid_pairs <- rownames(
    pair_table
  )[
    pair_table[, "Untreated"] == 1 &
      pair_table[, "Treated"] == 1
  ]
  
  paired_metadata <- paired_metadata %>%
    filter(
      biological_subject %in%
        valid_pairs
    ) %>%
    droplevels()
  
  n_pairs <- n_distinct(
    paired_metadata$biological_subject
  )
  
  if (
    n_pairs <
    CONFIG$min_paired_subjects
  ) {
    
    return(
      list(
        status = "Skipped",
        reason = paste0(
          "Only ",
          n_pairs,
          " exact complete pairs"
        ),
        result = NULL
      )
    )
  }
  
  paired_counts <- pb_counts[
    ,
    paired_metadata$pseudobulk_id,
    drop = FALSE
  ]
  
  design <- model.matrix(
    ~ biological_subject + condition,
    data = paired_metadata
  )
  
  dge <- DGEList(
    counts = paired_counts
  )
  
  keep <- filterByExpr(
    dge,
    design = design,
    min.count =
      CONFIG$min_count,
    min.total.count =
      CONFIG$min_total_count
  )
  
  if (
    sum(keep) <
    CONFIG$min_genes_after_filtering
  ) {
    
    return(
      list(
        status = "Skipped",
        reason = paste0(
          "Only ",
          sum(keep),
          " genes passed paired filtering"
        ),
        result = NULL
      )
    )
  }
  
  dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  dge <- calcNormFactors(
    dge,
    method = "TMM"
  )
  
  voom_object <- voom(
    dge,
    design,
    plot = FALSE
  )
  
  fit <- lmFit(
    voom_object,
    design
  )
  
  fit <- eBayes(
    fit,
    robust =
      CONFIG$robust_ebayes
  )
  
  treatment_coefficient <- grep(
    "^conditionTreated$",
    colnames(design),
    value = TRUE
  )
  
  if (
    length(
      treatment_coefficient
    ) != 1
  ) {
    
    return(
      list(
        status = "Failed",
        reason = paste0(
          "Could not find conditionTreated coefficient. Design columns: ",
          paste(
            colnames(design),
            collapse = ", "
          )
        ),
        result = NULL
      )
    )
  }
  
  result <- topTable(
    fit,
    coef =
      treatment_coefficient,
    number = Inf,
    sort.by = "P"
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    mutate(
      contrast =
        "Treated_vs_Untreated",
      numerator =
        "Treated",
      denominator =
        "Untreated",
      paired = TRUE,
      n_paired_subjects =
        n_pairs,
      celltype =
        celltype_value,
      celltype_clean =
        celltype_clean
    )
  
  list(
    status = "Completed",
    reason = NA_character_,
    result = result,
    voom = voom_object,
    fit = fit,
    design = design,
    paired_metadata =
      paired_metadata
  )
}

# =============================================================================
# SAVE ONE DEG RESULT AND RANKED FILE
# =============================================================================

save_deg_result <- function(
    result,
    celltype_clean
) {
  
  contrast_name <- unique(
    result$contrast
  )
  
  celltype_deg_dir <- file.path(
    DEG_DIR,
    celltype_clean
  )
  
  dir_create(
    celltype_deg_dir
  )
  
  deg_file <- file.path(
    celltype_deg_dir,
    paste0(
      "DEG_",
      celltype_clean,
      "_",
      contrast_name,
      ".csv"
    )
  )
  
  fwrite(
    result,
    deg_file
  )
  
  ranked <- result %>%
    transmute(
      gene,
      statistic = t,
      logFC,
      P.Value,
      adj.P.Val
    ) %>%
    filter(
      is.finite(statistic)
    ) %>%
    arrange(
      desc(statistic)
    )
  
  ranked_file <- file.path(
    celltype_deg_dir,
    paste0(
      "RANKED_",
      celltype_clean,
      "_",
      contrast_name,
      ".csv"
    )
  )
  
  fwrite(
    ranked,
    ranked_file
  )
  
  if (
    isTRUE(
      CONFIG$make_volcano_plots
    )
  ) {
    
    volcano <- make_volcano(
      result = result,
      contrast =
        contrast_name,
      celltype =
        unique(result$celltype)
    )
    
    save_plot_both(
      volcano,
      file.path(
        FIGURE_DIR,
        paste0(
          "Volcano_",
          celltype_clean,
          "_",
          contrast_name
        )
      ),
      width = 8,
      height = 6
    )
  }
  
  list(
    deg_file = deg_file,
    ranked_file = ranked_file
  )
}

# =============================================================================
# ANALYZE ONE CELL TYPE
# =============================================================================

analyse_celltype_deg <- function(
    celltype_value
) {
  
  celltype_clean <- clean_name(
    celltype_value
  )
  
  logmsg(
    "DEG: %s",
    celltype_value
  )
  
  ct_metadata <- metadata %>%
    filter(
      celltype ==
        celltype_value
    )
  
  ct_counts_summary <- ct_metadata %>%
    count(
      pseudobulk_id,
      dataset,
      biological_subject,
      subject_uid,
      condition,
      name = "n_cells"
    )
  
  valid_pseudobulks <- ct_counts_summary %>%
    filter(
      n_cells >=
        CONFIG$min_cells_per_pseudobulk
    ) %>%
    pull(
      pseudobulk_id
    )
  
  ct_metadata <- ct_metadata %>%
    filter(
      pseudobulk_id %in%
        valid_pseudobulks
    )
  
  if (nrow(ct_metadata) == 0) {
    
    return(
      list(
        celltype = celltype_value,
        status = "Skipped",
        reason =
          "No pseudobulks passed the minimum-cell threshold",
        contrast_status = NULL
      )
    )
  }
  
  ct_counts <- counts_matrix[
    ,
    ct_metadata$cell_barcode,
    drop = FALSE
  ]
  
  pb_counts <- collapse_sparse_columns(
    ct_counts,
    ct_metadata$pseudobulk_id
  )
  
  pb_metadata <- ct_metadata %>%
    distinct(
      pseudobulk_id,
      dataset,
      biological_subject,
      subject_uid,
      condition
    )
  
  pb_metadata <- pb_metadata[
    match(
      colnames(pb_counts),
      pb_metadata$pseudobulk_id
    ),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      pb_metadata$pseudobulk_id,
      colnames(pb_counts)
    )
  )
  
  pb_metadata <- pb_metadata %>%
    left_join(
      ct_counts_summary %>%
        select(
          pseudobulk_id,
          n_cells
        ),
      by = "pseudobulk_id"
    )
  
  celltype_pb_dir <- file.path(
    PB_DIR,
    celltype_clean
  )
  
  dir_create(
    celltype_pb_dir
  )
  
  write_csv(
    pb_metadata,
    file.path(
      celltype_pb_dir,
      paste0(
        celltype_clean,
        "_pseudobulk_metadata.csv"
      )
    )
  )
  
  if (
    isTRUE(
      CONFIG$save_pseudobulk_counts
    )
  ) {
    
    fwrite(
      as.data.frame(
        pb_counts
      ) %>%
        rownames_to_column(
          "gene"
        ),
      file.path(
        celltype_pb_dir,
        paste0(
          celltype_clean,
          "_raw_pseudobulk_counts.csv.gz"
        )
      )
    )
  }
  
  contrast_records <- list()
  saved_results <- list()
  
  # ---------------------------------------------------------------------------
  # Unpaired Untreated vs Healthy
  # ---------------------------------------------------------------------------
  
  untreated_healthy <-
    fit_unpaired_contrast(
      pb_counts =
        pb_counts,
      pb_metadata =
        pb_metadata,
      numerator =
        "Untreated",
      denominator =
        "Healthy",
      celltype_value =
        celltype_value,
      celltype_clean =
        celltype_clean
    )
  
  contrast_records[["Untreated_vs_Healthy"]] <- tibble(
    celltype =
      celltype_value,
    contrast =
      "Untreated_vs_Healthy",
    design =
      "Unpaired",
    status =
      untreated_healthy$status,
    reason =
      untreated_healthy$reason
  )
  
  if (
    identical(
      untreated_healthy$status,
      "Completed"
    )
  ) {
    
    save_deg_result(
      untreated_healthy$result,
      celltype_clean
    )
    
    saved_results[["Untreated_vs_Healthy"]] <- untreated_healthy$result
  }
  
  # ---------------------------------------------------------------------------
  # Unpaired Treated vs Healthy
  # ---------------------------------------------------------------------------
  
  treated_healthy <-
    fit_unpaired_contrast(
      pb_counts =
        pb_counts,
      pb_metadata =
        pb_metadata,
      numerator =
        "Treated",
      denominator =
        "Healthy",
      celltype_value =
        celltype_value,
      celltype_clean =
        celltype_clean
    )
  
  contrast_records[["Treated_vs_Healthy"]] <- tibble(
    celltype =
      celltype_value,
    contrast =
      "Treated_vs_Healthy",
    design =
      "Unpaired",
    status =
      treated_healthy$status,
    reason =
      treated_healthy$reason
  )
  
  if (
    identical(
      treated_healthy$status,
      "Completed"
    )
  ) {
    
    save_deg_result(
      treated_healthy$result,
      celltype_clean
    )
    
    saved_results[["Treated_vs_Healthy"]] <- treated_healthy$result
  }
  
  # ---------------------------------------------------------------------------
  # Explicit paired Treated vs Untreated
  # ---------------------------------------------------------------------------
  
  treated_untreated <-
    fit_paired_treatment_contrast(
      pb_counts =
        pb_counts,
      pb_metadata =
        pb_metadata,
      celltype_value =
        celltype_value,
      celltype_clean =
        celltype_clean
    )
  
  contrast_records[["Treated_vs_Untreated"]] <- tibble(
    celltype =
      celltype_value,
    contrast =
      "Treated_vs_Untreated",
    design =
      "Paired: ~ biological_subject + condition",
    status =
      treated_untreated$status,
    reason =
      treated_untreated$reason
  )
  
  if (
    identical(
      treated_untreated$status,
      "Completed"
    )
  ) {
    
    save_deg_result(
      treated_untreated$result,
      celltype_clean
    )
    
    saved_results[["Treated_vs_Untreated"]] <- treated_untreated$result
    
    write_csv(
      treated_untreated$paired_metadata,
      file.path(
        celltype_pb_dir,
        paste0(
          celltype_clean,
          "_paired_subject_metadata.csv"
        )
      )
    )
  }
  
  contrast_status <- bind_rows(
    contrast_records
  )
  
  write_csv(
    contrast_status,
    file.path(
      celltype_pb_dir,
      paste0(
        celltype_clean,
        "_contrast_status.csv"
      )
    )
  )
  
  if (
    isTRUE(
      CONFIG$save_model_objects
    )
  ) {
    
    saveRDS(
      list(
        Untreated_vs_Healthy =
          untreated_healthy,
        Treated_vs_Healthy =
          treated_healthy,
        Treated_vs_Untreated =
          treated_untreated
      ),
      file.path(
        celltype_pb_dir,
        paste0(
          celltype_clean,
          "_limma_models.rds"
        )
      )
    )
  }
  
  list(
    celltype = celltype_value,
    status = "Completed",
    reason = NA_character_,
    contrast_status =
      contrast_status,
    results =
      saved_results
  )
}

# =============================================================================
# RUN DEG ANALYSIS
# =============================================================================

celltypes <- metadata %>%
  distinct(
    celltype
  ) %>%
  arrange(
    celltype
  ) %>%
  pull(
    celltype
  )

logmsg(
  "Cell types discovered: %d",
  length(celltypes)
)

safe_analyse_deg <- safely(
  analyse_celltype_deg,
  otherwise = NULL,
  quiet = FALSE
)

deg_safe_results <- map(
  celltypes,
  safe_analyse_deg
)

names(deg_safe_results) <- celltypes

deg_analysis_results <- map(
  deg_safe_results,
  "result"
)

deg_unexpected_errors <- imap_dfr(
  deg_safe_results,
  function(x, celltype_name) {
    
    if (is.null(x$error)) {
      return(NULL)
    }
    
    tibble(
      celltype =
        celltype_name,
      error =
        conditionMessage(
          x$error
        )
    )
  }
)

if (
  nrow(
    deg_unexpected_errors
  ) > 0
) {
  
  write_csv(
    deg_unexpected_errors,
    file.path(
      SUMMARY_DIR,
      "unexpected_deg_errors.csv"
    )
  )
}

all_contrast_status <- deg_analysis_results %>%
  compact() %>%
  map(
    "contrast_status"
  ) %>%
  compact() %>%
  bind_rows()

write_csv(
  all_contrast_status,
  file.path(
    SUMMARY_DIR,
    "all_celltype_contrast_status.csv"
  )
)

print(
  all_contrast_status %>%
    count(
      contrast,
      design,
      status,
      name = "n_celltypes"
    ),
  n = Inf
)

# =============================================================================
# DEG COUNT SUMMARY
# =============================================================================

deg_files <- dir_ls(
  DEG_DIR,
  recurse = TRUE,
  regexp = "DEG_.*\\.csv$"
)

deg_count_summary <- map_dfr(
  deg_files,
  function(path) {
    
    x <- fread(
      path,
      data.table = FALSE
    )
    
    if (nrow(x) == 0) {
      return(NULL)
    }
    
    tibble(
      celltype =
        unique(x$celltype)[1],
      contrast =
        unique(x$contrast)[1],
      paired =
        unique(x$paired)[1],
      n_up = sum(
        x$adj.P.Val <=
          CONFIG$deg_fdr &
          x$logFC > 0,
        na.rm = TRUE
      ),
      n_down = sum(
        x$adj.P.Val <=
          CONFIG$deg_fdr &
          x$logFC < 0,
        na.rm = TRUE
      ),
      n_tested =
        nrow(x)
    ) %>%
      mutate(
        n_total =
          n_up + n_down
      )
  }
)

write_csv(
  deg_count_summary,
  file.path(
    SUMMARY_DIR,
    "DEG_count_summary.csv"
  )
)

# =============================================================================
# MSIGDB RETRIEVAL
# =============================================================================

get_msigdb <- function(
    collection,
    subcollection = NULL
) {
  
  available_arguments <- names(
    formals(
      msigdbr::msigdbr
    )
  )
  
  if (
    "collection" %in%
    available_arguments
  ) {
    
    arguments <- list(
      species =
        "Homo sapiens",
      collection =
        collection
    )
    
    if (!is.null(subcollection)) {
      arguments$subcollection <-
        subcollection
    }
    
  } else {
    
    arguments <- list(
      species =
        "Homo sapiens",
      category =
        collection
    )
    
    if (!is.null(subcollection)) {
      arguments$subcategory <-
        subcollection
    }
  }
  
  do.call(
    msigdbr::msigdbr,
    arguments
  )
}

msig_to_list <- function(msig_table) {
  
  set_column <- if (
    "gs_name" %in%
    colnames(msig_table)
  ) {
    "gs_name"
  } else {
    "gs_id"
  }
  
  gene_column <- if (
    "gene_symbol" %in%
    colnames(msig_table)
  ) {
    "gene_symbol"
  } else {
    "human_gene_symbol"
  }
  
  msig_table %>%
    filter(
      !is.na(
        .data[[set_column]]
      ),
      !is.na(
        .data[[gene_column]]
      )
    ) %>%
    distinct(
      .data[[set_column]],
      .data[[gene_column]]
    ) %>%
    split(
      x = .[[gene_column]],
      f = .[[set_column]]
    )
}

logmsg(
  "Loading Hallmark pathways"
)

hallmark <- get_msigdb(
  "H"
)

logmsg(
  "Loading Reactome pathways"
)

reactome <- get_msigdb(
  "C2",
  "CP:REACTOME"
)

logmsg(
  "Loading GO Biological Process pathways"
)

go_bp <- get_msigdb(
  "C5",
  "GO:BP"
)

PATHWAYS <- list(
  Hallmark =
    msig_to_list(hallmark),
  Reactome =
    msig_to_list(reactome),
  GO_BP =
    msig_to_list(go_bp)
)

# =============================================================================
# PREPARE RANKED VECTOR
# =============================================================================

prepare_ranked_vector <- function(
    ranked_file
) {
  
  ranked_table <- fread(
    ranked_file,
    data.table = FALSE
  ) %>%
    transmute(
      gene =
        as.character(gene),
      statistic =
        as.numeric(statistic)
    ) %>%
    filter(
      !is.na(gene),
      gene != "",
      is.finite(statistic)
    ) %>%
    group_by(
      gene
    ) %>%
    slice_max(
      order_by =
        abs(statistic),
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    arrange(
      desc(statistic)
    )
  
  stats <- ranked_table$statistic
  
  names(stats) <-
    ranked_table$gene
  
  stats
}

# =============================================================================
# RUN ONE FGSEA ANALYSIS
# =============================================================================

run_fgsea_one <- function(
    ranked_file,
    database
) {
  
  stats <- prepare_ranked_vector(
    ranked_file
  )
  
  pathways <- map(
    PATHWAYS[[database]],
    intersect,
    names(stats)
  )
  
  pathways <- pathways[
    lengths(pathways) >=
      CONFIG$gsea_min_size &
      lengths(pathways) <=
      CONFIG$gsea_max_size
  ]
  
  if (
    length(pathways) == 0
  ) {
    return(NULL)
  }
  
  result <- fgseaMultilevel(
    pathways =
      pathways,
    stats =
      stats,
    minSize =
      CONFIG$gsea_min_size,
    maxSize =
      CONFIG$gsea_max_size,
    eps =
      CONFIG$gsea_eps,
    scoreType =
      "std"
  )
  
  as.data.frame(result)
}

# =============================================================================
# PATHWAY PLOT
# =============================================================================

make_pathway_plot <- function(
    result,
    celltype,
    contrast,
    database
) {
  
  plot_data <- result %>%
    filter(
      !is.na(padj),
      padj <=
        CONFIG$pathway_fdr
    )
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  positive <- plot_data %>%
    filter(NES > 0) %>%
    arrange(
      padj,
      desc(NES)
    ) %>%
    slice_head(
      n =
        CONFIG$top_n_pathways
    )
  
  negative <- plot_data %>%
    filter(NES < 0) %>%
    arrange(
      padj,
      NES
    ) %>%
    slice_head(
      n =
        CONFIG$top_n_pathways
    )
  
  plot_data <- bind_rows(
    positive,
    negative
  ) %>%
    distinct(
      pathway,
      .keep_all = TRUE
    ) %>%
    mutate(
      pathway_label =
        pathway %>%
        str_remove(
          "^HALLMARK_"
        ) %>%
        str_remove(
          "^REACTOME_"
        ) %>%
        str_remove(
          "^GOBP_"
        ) %>%
        str_replace_all(
          "_",
          " "
        ) %>%
        str_to_title() %>%
        str_wrap(
          width = 45
        ),
      
      pathway_label =
        fct_reorder(
          pathway_label,
          NES
        ),
      
      direction =
        if_else(
          NES > 0,
          "Positive NES",
          "Negative NES"
        )
    )
  
  ggplot(
    plot_data,
    aes(
      x = NES,
      y = pathway_label,
      fill = direction
    )
  ) +
    geom_col(
      width = 0.75
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title = paste0(
        celltype,
        ": ",
        database
      ),
      subtitle =
        contrast,
      x =
        "Normalized enrichment score",
      y = NULL,
      fill = NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      plot.title =
        element_text(
          face = "bold"
        ),
      panel.grid.minor =
        element_blank(),
      legend.position =
        "bottom"
    )
}

# =============================================================================
# DISCOVER AND RUN ALL RANKED FILES
# =============================================================================

ranked_files <- dir_ls(
  DEG_DIR,
  recurse = TRUE,
  regexp = "RANKED_.*\\.csv$"
)

if (
  length(ranked_files) == 0
) {
  
  stop(
    "No ranked files were generated."
  )
}

parse_ranked_file <- function(path) {
  
  filename <- basename(path)
  
  contrast <- case_when(
    
    str_detect(
      filename,
      "Untreated_vs_Healthy"
    ) ~
      "Untreated_vs_Healthy",
    
    str_detect(
      filename,
      "Treated_vs_Healthy"
    ) ~
      "Treated_vs_Healthy",
    
    str_detect(
      filename,
      "Treated_vs_Untreated"
    ) ~
      "Treated_vs_Untreated",
    
    TRUE ~
      NA_character_
  )
  
  tibble(
    ranked_file = path,
    celltype_clean =
      basename(dirname(path)),
    contrast = contrast
  )
}

ranked_manifest <- map_dfr(
  ranked_files,
  parse_ranked_file
) %>%
  filter(
    !is.na(contrast)
  )

gsea_grid <- crossing(
  ranked_manifest,
  database =
    names(PATHWAYS)
)

safe_fgsea <- safely(
  function(
    ranked_file,
    celltype_clean,
    contrast,
    database
  ) {
    
    logmsg(
      "GSEA: %s | %s | %s",
      celltype_clean,
      contrast,
      database
    )
    
    result <- run_fgsea_one(
      ranked_file,
      database
    )
    
    if (
      is.null(result) ||
      nrow(result) == 0
    ) {
      
      return(NULL)
    }
    
    result <- result %>%
      mutate(
        leadingEdge =
          map_chr(
            leadingEdge,
            ~ paste(
              .x,
              collapse = ";"
            )
          ),
        
        celltype =
          celltype_clean,
        
        contrast =
          contrast,
        
        database =
          database,
        
        significant =
          !is.na(padj) &
          padj <=
          CONFIG$pathway_fdr
      ) %>%
      arrange(
        padj,
        desc(abs(NES))
      )
    
    output_dir <- file.path(
      GSEA_DIR,
      celltype_clean
    )
    
    dir_create(
      output_dir
    )
    
    output_file <- file.path(
      output_dir,
      paste0(
        celltype_clean,
        "_",
        contrast,
        "_",
        database,
        "_all.csv"
      )
    )
    
    significant_file <- file.path(
      output_dir,
      paste0(
        celltype_clean,
        "_",
        contrast,
        "_",
        database,
        "_significant_FDR0p10.csv"
      )
    )
    
    fwrite(
      result,
      output_file
    )
    
    fwrite(
      result %>%
        filter(significant),
      significant_file
    )
    
    if (
      isTRUE(
        CONFIG$make_pathway_plots
      )
    ) {
      
      pathway_plot <-
        make_pathway_plot(
          result =
            result,
          celltype =
            celltype_clean,
          contrast =
            contrast,
          database =
            database
        )
      
      if (!is.null(pathway_plot)) {
        
        save_plot_both(
          pathway_plot,
          file.path(
            FIGURE_DIR,
            paste0(
              "TopPathways_",
              celltype_clean,
              "_",
              contrast,
              "_",
              database
            )
          ),
          width = 12,
          height = 10
        )
      }
    }
    
    result
  },
  otherwise = NULL,
  quiet = FALSE
)

gsea_safe_results <- pmap(
  gsea_grid,
  safe_fgsea
)

all_gsea_results <- map(
  gsea_safe_results,
  "result"
) %>%
  compact() %>%
  compact() %>%
  bind_rows()

if (
  nrow(all_gsea_results) == 0
) {
  
  stop(
    "No GSEA results were produced."
  )
}

write_csv(
  all_gsea_results,
  file.path(
    SUMMARY_DIR,
    "all_GSEA_results.csv"
  )
)

write_csv(
  all_gsea_results %>%
    filter(
      padj <=
        CONFIG$pathway_fdr
    ),
  file.path(
    SUMMARY_DIR,
    "all_significant_GSEA_results_FDR0p10.csv"
  )
)

# =============================================================================
# PATHWAY COUNT SUMMARY
# =============================================================================

pathway_count_summary <- all_gsea_results %>%
  group_by(
    celltype,
    contrast,
    database
  ) %>%
  summarise(
    n_tested = n(),
    
    n_significant = sum(
      padj <=
        CONFIG$pathway_fdr,
      na.rm = TRUE
    ),
    
    n_positive = sum(
      padj <=
        CONFIG$pathway_fdr &
        NES > 0,
      na.rm = TRUE
    ),
    
    n_negative = sum(
      padj <=
        CONFIG$pathway_fdr &
        NES < 0,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

write_csv(
  pathway_count_summary,
  file.path(
    SUMMARY_DIR,
    "pathway_count_summary.csv"
  )
)

# =============================================================================
# JOIN THREE CONTRASTS
# =============================================================================

prepare_for_join <- function(
    contrast_name,
    suffix
) {
  
  all_gsea_results %>%
    filter(
      contrast ==
        contrast_name
    ) %>%
    select(
      celltype,
      database,
      pathway,
      NES,
      padj,
      pval,
      leadingEdge
    ) %>%
    rename(
      !!paste0(
        "NES_",
        suffix
      ) := NES,
      
      !!paste0(
        "padj_",
        suffix
      ) := padj,
      
      !!paste0(
        "pval_",
        suffix
      ) := pval,
      
      !!paste0(
        "leadingEdge_",
        suffix
      ) := leadingEdge
    )
}

UH <- prepare_for_join(
  "Untreated_vs_Healthy",
  "UH"
)

TH <- prepare_for_join(
  "Treated_vs_Healthy",
  "TH"
)

TU <- prepare_for_join(
  "Treated_vs_Untreated",
  "TU"
)

joined_pathways <- full_join(
  UH,
  TH,
  by = c(
    "celltype",
    "database",
    "pathway"
  )
) %>%
  full_join(
    TU,
    by = c(
      "celltype",
      "database",
      "pathway"
    )
  )

# =============================================================================
# THREE-CONTRAST PATHWAY CLASSIFICATION
#
# UH = Untreated vs Healthy
# TH = Treated vs Healthy
# TU = Treated vs Untreated, explicitly paired
# =============================================================================

FDR <- CONFIG$pathway_fdr

classified_pathways <- joined_pathways %>%
  mutate(
    
    disease_sig =
      !is.na(padj_UH) &
      padj_UH <= FDR,
    
    treated_vs_healthy_sig =
      !is.na(padj_TH) &
      padj_TH <= FDR,
    
    paired_treatment_sig =
      !is.na(padj_TU) &
      padj_TU <= FDR,
    
    moves_toward_healthy =
      disease_sig &
      paired_treatment_sig &
      !is.na(NES_UH) &
      !is.na(NES_TU) &
      sign(NES_UH) ==
      -sign(NES_TU),
    
    moves_away_from_healthy =
      disease_sig &
      paired_treatment_sig &
      !is.na(NES_UH) &
      !is.na(NES_TU) &
      sign(NES_UH) ==
      sign(NES_TU),
    
    residual_same_direction =
      disease_sig &
      treated_vs_healthy_sig &
      !is.na(NES_UH) &
      !is.na(NES_TH) &
      sign(NES_UH) ==
      sign(NES_TH),
    
    residual_opposite_direction =
      disease_sig &
      treated_vs_healthy_sig &
      !is.na(NES_UH) &
      !is.na(NES_TH) &
      sign(NES_UH) ==
      -sign(NES_TH),
    
    residual_ratio = if_else(
      !is.na(NES_UH) &
        abs(NES_UH) > 0 &
        !is.na(NES_TH),
      
      abs(NES_TH) /
        abs(NES_UH),
      
      NA_real_
    ),
    
    pathway_category = case_when(
      
      disease_sig &
        moves_toward_healthy &
        !treated_vs_healthy_sig ~
        "Recovered toward healthy",
      
      disease_sig &
        moves_toward_healthy &
        residual_same_direction &
        !is.na(residual_ratio) &
        residual_ratio < 1 ~
        "Partially recovered",
      
      disease_sig &
        moves_toward_healthy &
        residual_opposite_direction ~
        "Over-corrected beyond healthy",
      
      disease_sig &
        moves_away_from_healthy ~
        "Moved away from healthy",
      
      disease_sig &
        residual_same_direction &
        !moves_toward_healthy ~
        "Persistent disease pathway",
      
      disease_sig &
        !moves_toward_healthy &
        !moves_away_from_healthy ~
        "No clear paired treatment response",
      
      !disease_sig &
        paired_treatment_sig ~
        "Treatment-induced pathway",
      
      TRUE ~
        "Not classified"
    ),
    
    disease_direction = case_when(
      NES_UH > 0 ~
        "Higher in untreated than healthy",
      NES_UH < 0 ~
        "Lower in untreated than healthy",
      TRUE ~
        NA_character_
    ),
    
    paired_treatment_direction =
      case_when(
        NES_TU > 0 ~
          "Higher after treatment",
        NES_TU < 0 ~
          "Lower after treatment",
        TRUE ~
          NA_character_
      )
  )

write_csv(
  joined_pathways,
  file.path(
    DECISION_DIR,
    "joined_three_contrasts.csv"
  )
)

write_csv(
  classified_pathways,
  file.path(
    DECISION_DIR,
    "three_contrast_pathway_classification.csv"
  )
)

write_csv(
  classified_pathways %>%
    filter(disease_sig),
  file.path(
    DECISION_DIR,
    "disease_pathway_recovery_classification.csv"
  )
)

# =============================================================================
# CATEGORY-SPECIFIC OUTPUTS
# =============================================================================

categories <- c(
  "Recovered toward healthy",
  "Partially recovered",
  "Over-corrected beyond healthy",
  "Persistent disease pathway",
  "Moved away from healthy",
  "No clear paired treatment response",
  "Treatment-induced pathway"
)

walk(
  categories,
  function(category_name) {
    
    filename <- clean_name(
      tolower(category_name)
    )
    
    write_csv(
      classified_pathways %>%
        filter(
          pathway_category ==
            category_name
        ),
      file.path(
        DECISION_DIR,
        paste0(
          filename,
          ".csv"
        )
      )
    )
  }
)

# =============================================================================
# CATEGORY COUNT SUMMARY
# =============================================================================

category_counts <- classified_pathways %>%
  filter(
    pathway_category !=
      "Not classified"
  ) %>%
  count(
    celltype,
    database,
    pathway_category,
    name = "n_pathways"
  )

write_csv(
  category_counts,
  file.path(
    DECISION_DIR,
    "pathway_category_counts.csv"
  )
)

grand_category_counts <- classified_pathways %>%
  filter(
    pathway_category !=
      "Not classified"
  ) %>%
  count(
    pathway_category,
    name = "n_pathways"
  ) %>%
  arrange(
    desc(n_pathways)
  )

write_csv(
  grand_category_counts,
  file.path(
    DECISION_DIR,
    "grand_pathway_category_counts.csv"
  )
)

print(
  grand_category_counts,
  n = Inf
)

# =============================================================================
# TOP RECOVERED AND PERSISTENT PATHWAYS
# =============================================================================

top_recovered <- classified_pathways %>%
  filter(
    pathway_category %in%
      c(
        "Recovered toward healthy",
        "Partially recovered"
      )
  ) %>%
  arrange(
    celltype,
    database,
    padj_UH,
    padj_TU
  ) %>%
  group_by(
    celltype,
    database
  ) %>%
  slice_head(
    n = 20
  ) %>%
  ungroup()

write_csv(
  top_recovered,
  file.path(
    DECISION_DIR,
    "top_recovered_pathways.csv"
  )
)

top_persistent <- classified_pathways %>%
  filter(
    pathway_category ==
      "Persistent disease pathway"
  ) %>%
  arrange(
    celltype,
    database,
    padj_UH,
    padj_TH
  ) %>%
  group_by(
    celltype,
    database
  ) %>%
  slice_head(
    n = 20
  ) %>%
  ungroup()

write_csv(
  top_persistent,
  file.path(
    DECISION_DIR,
    "top_persistent_pathways.csv"
  )
)

# =============================================================================
# FINAL COMPLETION SUMMARY
# =============================================================================

completion_summary <- tibble(
  metric = c(
    "Cell types discovered",
    "Completed DEG contrasts",
    "Skipped DEG contrasts",
    "DEG files generated",
    "Ranked files generated",
    "Total GSEA results",
    "Significant GSEA results FDR 0.10",
    "Disease pathways classified",
    "Recovered pathways",
    "Partially recovered pathways",
    "Persistent disease pathways"
  ),
  
  value = c(
    length(celltypes),
    
    sum(
      all_contrast_status$status ==
        "Completed",
      na.rm = TRUE
    ),
    
    sum(
      all_contrast_status$status ==
        "Skipped",
      na.rm = TRUE
    ),
    
    length(
      dir_ls(
        DEG_DIR,
        recurse = TRUE,
        regexp = "DEG_.*\\.csv$"
      )
    ),
    
    length(
      ranked_files
    ),
    
    nrow(
      all_gsea_results
    ),
    
    sum(
      all_gsea_results$padj <=
        CONFIG$pathway_fdr,
      na.rm = TRUE
    ),
    
    sum(
      classified_pathways$disease_sig,
      na.rm = TRUE
    ),
    
    sum(
      classified_pathways$pathway_category ==
        "Recovered toward healthy",
      na.rm = TRUE
    ),
    
    sum(
      classified_pathways$pathway_category ==
        "Partially recovered",
      na.rm = TRUE
    ),
    
    sum(
      classified_pathways$pathway_category ==
        "Persistent disease pathway",
      na.rm = TRUE
    )
  )
)

write_csv(
  completion_summary,
  file.path(
    SUMMARY_DIR,
    "pipeline_completion_summary.csv"
  )
)

print(
  completion_summary,
  n = Inf
)

# =============================================================================
# SESSION INFO
# =============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    CONFIG$out_dir,
    "sessionInfo.txt"
  )
)

# =============================================================================
# FINAL REPORT
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("COMBINED LIMMA-VOOM + FGSEA PIPELINE COMPLETE\n")
cat("====================================================================\n")

cat("\nDEG contrasts:\n")
cat("  Untreated_vs_Healthy: unpaired\n")
cat("  Treated_vs_Healthy: unpaired\n")
cat(
  "  Treated_vs_Untreated: explicitly paired using ",
  "~ biological_subject + condition\n",
  sep = ""
)

cat("\nBiological replicate definitions:\n")
cat("  Yu2021: donorID_unified\n")
cat("  Internal: Subject_fixed\n")

cat("\nPathway collections:\n")
cat("  Hallmark\n")
cat("  Reactome\n")
cat("  GO Biological Process\n")

cat("\nMain DEG directory:\n")
cat(
  "  ",
  normalizePath(
    DEG_DIR,
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat("\nMain GSEA directory:\n")
cat(
  "  ",
  normalizePath(
    GSEA_DIR,
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat("\nMain pathway classification:\n")
cat(
  "  ",
  file.path(
    DECISION_DIR,
    "three_contrast_pathway_classification.csv"
  ),
  "\n",
  sep = ""
)

cat("\nPaired treatment result:\n")
cat(
  "  Treated_vs_Untreated uses only internal subjects with both ",
  "Untreated and Treated pseudobulks for the cell type.\n",
  sep = ""
)

cat("\nCross-study limitation:\n")
cat(
  "  Healthy occurs only in Yu2021. Untreated_vs_Healthy and ",
  "Treated_vs_Healthy remain cross-study comparisons.\n",
  sep = ""
)

cat("====================================================================\n")