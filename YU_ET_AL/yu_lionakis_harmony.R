# =============================================================================
# YU2021 + INTERNAL DUODENUM:
# QC METADATA, INDEPENDENT UMAPS, AND DATASET-LEVEL HARMONY INTEGRATION
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(readr)
  library(Matrix)
})

set.seed(1234)

# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG <- list(
  
  # ---------------------------------------------------------------------------
  # Input files
  # ---------------------------------------------------------------------------
  
  yu_rds = paste0(
    "YU2021_CONTROL_DUODENUM_REFERENCE/",
    "Yu2021_control_duodenum_live_biopsy.rds"
  ),
  
  internal_rds = paste0(
    "/vf/users/cheungf/ACEPD/",
    "LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/",
    "object_final_annotated_panGI_withDay.rds"
  ),
  
  # ---------------------------------------------------------------------------
  # Output
  # ---------------------------------------------------------------------------
  
  out_dir = "YU2021_INTERNAL_HARMONY",
  
  # ---------------------------------------------------------------------------
  # Assay and processing parameters
  # ---------------------------------------------------------------------------
  
  assay = "RNA",
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  n_variable_features = 3000,
  n_integration_features = 3000,
  n_pcs = 50,
  dims_use = 1:30,
  clustering_resolution = 0.5,
  
  # ---------------------------------------------------------------------------
  # Internal conditions used in the combined integration
  #
  # Spike is retained in the independent internal UMAP but excluded from
  # the primary Yu healthy versus treated/untreated integration.
  # ---------------------------------------------------------------------------
  
  internal_conditions_integrate = c(
    "Untreated",
    "Treated"
  ),
  
  # ---------------------------------------------------------------------------
  # Harmony
  #
  # Correct only the study/dataset source.
  # Do not correct Subject_fixed, donorID_unified, Treatment, or sampleID.
  # ---------------------------------------------------------------------------
  
  harmony_group = "dataset",
  
  # ---------------------------------------------------------------------------
  # Plot settings
  # ---------------------------------------------------------------------------
  
  raster = TRUE
)

dir.create(
  CONFIG$out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

comparison_dir <- file.path(
  CONFIG$out_dir,
  "01_QC_AND_COMPOSITION"
)

yu_dir <- file.path(
  CONFIG$out_dir,
  "02_YU2021_INDEPENDENT"
)

internal_dir <- file.path(
  CONFIG$out_dir,
  "03_INTERNAL_INDEPENDENT"
)

integration_dir <- file.path(
  CONFIG$out_dir,
  "04_HARMONY_INTEGRATION"
)

dir.create(
  comparison_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  yu_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  internal_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  integration_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# HELPER FUNCTIONS
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

theme_analysis <- function() {
  
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      plot.subtitle = element_text(
        size = 9
      ),
      axis.title = element_text(
        face = "bold"
      ),
      strip.text = element_text(
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(
        8,
        12,
        8,
        8
      )
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
    bg = "white"
  )
  
  ggsave(
    filename = paste0(
      filename_prefix,
      ".pdf"
    ),
    plot = plot_object,
    width = width,
    height = height,
    bg = "white"
  )
}

clean_character <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x[
    is.na(x) |
      x == ""
  ] <- NA_character_
  
  x
}

get_counts_matrix <- function(object) {
  
  tryCatch(
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
}

# =============================================================================
# ENSURE STANDARD QC METRICS
# =============================================================================

ensure_qc_metrics <- function(
    object,
    object_name
) {
  
  DefaultAssay(object) <- CONFIG$assay
  
  metadata_columns <- colnames(
    object@meta.data
  )
  
  need_counts <- (
    !"nCount_RNA" %in% metadata_columns
  ) ||
    all(
      is.na(
        object$nCount_RNA
      )
    )
  
  need_features <- (
    !"nFeature_RNA" %in% metadata_columns
  ) ||
    all(
      is.na(
        object$nFeature_RNA
      )
    )
  
  need_percent_mt <- (
    !"percent.mt" %in% metadata_columns
  ) ||
    all(
      is.na(
        object$percent.mt
      )
    )
  
  need_percent_ribo <- (
    !"percent.ribo" %in% metadata_columns
  ) ||
    all(
      is.na(
        object$percent.ribo
      )
    )
  
  need_percent_hb <- (
    !"percent.hb" %in% metadata_columns
  ) ||
    all(
      is.na(
        object$percent.hb
      )
    )
  
  if (need_counts || need_features) {
    
    logmsg(
      "%s: reading count matrix for QC calculations",
      object_name
    )
    
    counts_matrix <- get_counts_matrix(
      object
    )
    
    if (need_counts) {
      
      object$nCount_RNA <- Matrix::colSums(
        counts_matrix
      )
    }
    
    if (need_features) {
      
      object$nFeature_RNA <- Matrix::colSums(
        counts_matrix > 0
      )
    }
  }
  
  if (need_percent_mt) {
    
    object$percent.mt <- PercentageFeatureSet(
      object,
      assay = CONFIG$assay,
      pattern = "^MT-"
    )
  }
  
  if (need_percent_ribo) {
    
    object$percent.ribo <- PercentageFeatureSet(
      object,
      assay = CONFIG$assay,
      pattern = "^RP[SL][0-9]"
    )
  }
  
  if (need_percent_hb) {
    
    object$percent.hb <- PercentageFeatureSet(
      object,
      assay = CONFIG$assay,
      pattern = "^HB[ABDEGQZ][0-9]"
    )
  }
  
  qc_columns <- c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.ribo",
    "percent.hb"
  )
  
  for (column_name in qc_columns) {
    
    object@meta.data[[column_name]] <- as.numeric(
      object@meta.data[[column_name]]
    )
  }
  
  object
}

# =============================================================================
# LOAD DATA
# =============================================================================

logmsg(
  "STEP 1: Loading Yu2021 object"
)

if (!file.exists(CONFIG$yu_rds)) {
  
  stop(
    "Yu2021 RDS not found: ",
    CONFIG$yu_rds
  )
}

yu <- readRDS(
  CONFIG$yu_rds
)

logmsg(
  "  Yu2021: %s cells and %s genes",
  comma(ncol(yu)),
  comma(nrow(yu))
)

logmsg(
  "STEP 2: Loading internal object"
)

if (!file.exists(CONFIG$internal_rds)) {
  
  stop(
    "Internal RDS not found: ",
    CONFIG$internal_rds
  )
}

obj <- readRDS(
  CONFIG$internal_rds
)

logmsg(
  "  Internal: %s cells and %s genes",
  comma(ncol(obj)),
  comma(nrow(obj))
)

if (!inherits(yu, "Seurat")) {
  
  stop(
    "Yu input is not a Seurat object."
  )
}

if (!inherits(obj, "Seurat")) {
  
  stop(
    "Internal input is not a Seurat object."
  )
}

DefaultAssay(yu) <- CONFIG$assay
DefaultAssay(obj) <- CONFIG$assay

# =============================================================================
# VERIFY ORIGINAL METADATA
# =============================================================================

required_yu_columns <- c(
  "sampleID",
  "donorID_unified",
  "level_1_annot",
  "level_2_annot",
  "level_3_annot"
)

required_internal_columns <- c(
  "orig.ident",
  "Subject_fixed",
  "Treatment",
  "PanGI_L1",
  "PanGI_L2",
  "PanGI_L3"
)

missing_yu <- setdiff(
  required_yu_columns,
  colnames(yu@meta.data)
)

missing_internal <- setdiff(
  required_internal_columns,
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
# QC METRICS
# =============================================================================

logmsg(
  "STEP 3: Verifying QC metrics"
)

yu <- ensure_qc_metrics(
  yu,
  "Yu2021"
)

obj <- ensure_qc_metrics(
  obj,
  "Internal"
)

# =============================================================================
# STANDARDIZED METADATA
# =============================================================================

logmsg(
  "STEP 4: Creating standardized technical and biological metadata"
)

# -----------------------------------------------------------------------------
# Yu2021
# -----------------------------------------------------------------------------

yu$dataset <- "Yu2021"

# Technical sample identifier
yu$technical_sample <- clean_character(
  yu$sampleID
)

# Biological subject identifier
yu$biological_subject <- clean_character(
  yu$donorID_unified
)

yu$condition_common <- "Healthy"

yu$celltype_L1_common <- clean_character(
  yu$level_1_annot
)

yu$celltype_L2_common <- clean_character(
  yu$level_2_annot
)

yu$celltype_L3_common <- clean_character(
  yu$level_3_annot
)

# -----------------------------------------------------------------------------
# Internal
# -----------------------------------------------------------------------------

obj$dataset <- "Internal"

# Technical pooled GEX/library identifier
obj$technical_sample <- clean_character(
  obj$orig.ident
)

# Biological SNP-assigned subject identifier
obj$biological_subject <- clean_character(
  obj$Subject_fixed
)

obj$condition_common <- clean_character(
  obj$Treatment
)

obj$celltype_L1_common <- clean_character(
  obj$PanGI_L1
)

obj$celltype_L2_common <- clean_character(
  obj$PanGI_L2
)

obj$celltype_L3_common <- clean_character(
  obj$PanGI_L3
)

# =============================================================================
# MAKE CELL NAMES UNIQUE
# =============================================================================

if (!all(
  grepl(
    "^YU2021_",
    colnames(yu)
  )
)) {
  
  yu <- RenameCells(
    yu,
    add.cell.id = "YU2021"
  )
}

if (!all(
  grepl(
    "^INTERNAL_",
    colnames(obj)
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
    "Cell names overlap after prefixing."
  )
}

# =============================================================================
# SAMPLE AND SUBJECT MAPPING TABLES
# =============================================================================

yu_sample_subject <- yu@meta.data %>%
  distinct(
    technical_sample,
    biological_subject,
    condition_common
  ) %>%
  arrange(
    technical_sample
  )

internal_sample_subject <- obj@meta.data %>%
  distinct(
    technical_sample,
    biological_subject,
    condition_common
  ) %>%
  arrange(
    technical_sample,
    biological_subject
  )

write_csv(
  yu_sample_subject,
  file.path(
    comparison_dir,
    "Yu2021_sample_subject_mapping.csv"
  )
)

write_csv(
  internal_sample_subject,
  file.path(
    comparison_dir,
    "Internal_pool_subject_condition_mapping.csv"
  )
)

cat("\nYu2021 sample-to-donor mapping:\n")

print(
  table(
    yu$technical_sample,
    yu$biological_subject
  )
)

cat("\nInternal GEX pool-to-subject mapping:\n")

print(
  table(
    obj$technical_sample,
    obj$biological_subject
  )
)

cat("\nInternal subject-to-treatment mapping:\n")

print(
  table(
    obj$biological_subject,
    obj$condition_common
  )
)

# =============================================================================
# INDEPENDENT PROCESSING FUNCTION
# =============================================================================

process_independently <- function(
    object,
    object_name,
    output_directory
) {
  
  logmsg(
    "%s: beginning independent PCA and UMAP",
    object_name
  )
  
  DefaultAssay(object) <- CONFIG$assay
  
  # Remove standard reductions that will be recalculated.
  # Yu scVI and scANVI reductions are retained.
  reductions_remove <- intersect(
    c(
      "pca",
      "umap",
      "umap.pca",
      "harmony",
      "umap.harmony",
      "umap.preharmony"
    ),
    Reductions(object)
  )
  
  for (
    reduction_name in reductions_remove
  ) {
    
    object[[reduction_name]] <- NULL
  }
  
  existing_graphs <- names(
    object@graphs
  )
  
  if (length(existing_graphs) > 0) {
    
    for (
      graph_name in existing_graphs
    ) {
      
      object[[graph_name]] <- NULL
    }
  }
  
  object <- NormalizeData(
    object,
    assay = CONFIG$assay,
    normalization.method =
      CONFIG$normalization_method,
    scale.factor =
      CONFIG$scale_factor,
    verbose = FALSE
  )
  
  object <- FindVariableFeatures(
    object,
    assay = CONFIG$assay,
    selection.method = "vst",
    nfeatures =
      CONFIG$n_variable_features,
    verbose = FALSE
  )
  
  object <- ScaleData(
    object,
    assay = CONFIG$assay,
    features = VariableFeatures(object),
    verbose = FALSE
  )
  
  maximum_pcs <- min(
    CONFIG$n_pcs,
    length(
      VariableFeatures(object)
    ) - 1,
    ncol(object) - 1
  )
  
  object <- RunPCA(
    object,
    assay = CONFIG$assay,
    features = VariableFeatures(object),
    npcs = maximum_pcs,
    reduction.name = "pca",
    reduction.key = "PC_",
    seed.use = 1234,
    verbose = FALSE
  )
  
  dims_use <- CONFIG$dims_use[
    CONFIG$dims_use <= maximum_pcs
  ]
  
  object <- FindNeighbors(
    object,
    reduction = "pca",
    dims = dims_use,
    graph.name = c(
      "RNA_pca_nn",
      "RNA_pca_snn"
    ),
    verbose = FALSE
  )
  
  object <- FindClusters(
    object,
    graph.name = "RNA_pca_snn",
    resolution =
      CONFIG$clustering_resolution,
    cluster.name = "independent_cluster",
    random.seed = 1234,
    verbose = FALSE
  )
  
  object <- RunUMAP(
    object,
    reduction = "pca",
    dims = dims_use,
    reduction.name = "umap.pca",
    reduction.key = "pUMAP_",
    seed.use = 1234,
    verbose = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # UMAP by technical sample
  # ---------------------------------------------------------------------------
  
  p_technical <- DimPlot(
    object,
    reduction = "umap.pca",
    group.by = "technical_sample",
    raster = CONFIG$raster
  ) +
    ggtitle(
      paste0(
        object_name,
        ": PCA UMAP by technical sample"
      )
    )
  
  # ---------------------------------------------------------------------------
  # UMAP by biological subject
  # ---------------------------------------------------------------------------
  
  p_subject <- DimPlot(
    object,
    reduction = "umap.pca",
    group.by = "biological_subject",
    raster = CONFIG$raster
  ) +
    ggtitle(
      paste0(
        object_name,
        ": PCA UMAP by biological subject"
      )
    )
  
  # ---------------------------------------------------------------------------
  # UMAP by biological condition
  # ---------------------------------------------------------------------------
  
  p_condition <- DimPlot(
    object,
    reduction = "umap.pca",
    group.by = "condition_common",
    raster = CONFIG$raster
  ) +
    ggtitle(
      paste0(
        object_name,
        ": PCA UMAP by condition"
      )
    )
  
  # ---------------------------------------------------------------------------
  # UMAP by PanGI L2
  # ---------------------------------------------------------------------------
  
  p_celltype <- DimPlot(
    object,
    reduction = "umap.pca",
    group.by = "celltype_L2_common",
    label = TRUE,
    repel = TRUE,
    raster = CONFIG$raster
  ) +
    ggtitle(
      paste0(
        object_name,
        ": PCA UMAP by PanGI L2"
      )
    ) +
    NoLegend()
  
  # ---------------------------------------------------------------------------
  # Combined diagnostic figure
  # ---------------------------------------------------------------------------
  
  p_combined <- (
    p_technical |
      p_subject
  ) /
    (
      p_condition |
        p_celltype
    ) +
    plot_layout(
      guides = "collect"
    ) +
    plot_annotation(
      title = paste0(
        object_name,
        ": independent technical and biological diagnostics"
      )
    ) &
    theme(
      legend.position = "bottom"
    )
  
  save_plot_both(
    p_technical,
    file.path(
      output_directory,
      paste0(
        object_name,
        "_UMAP_technical_sample"
      )
    ),
    width = 11,
    height = 8
  )
  
  save_plot_both(
    p_subject,
    file.path(
      output_directory,
      paste0(
        object_name,
        "_UMAP_biological_subject"
      )
    ),
    width = 11,
    height = 8
  )
  
  save_plot_both(
    p_condition,
    file.path(
      output_directory,
      paste0(
        object_name,
        "_UMAP_condition"
      )
    ),
    width = 10,
    height = 8
  )
  
  save_plot_both(
    p_celltype,
    file.path(
      output_directory,
      paste0(
        object_name,
        "_UMAP_PanGI_L2"
      )
    ),
    width = 11,
    height = 8
  )
  
  save_plot_both(
    p_combined,
    file.path(
      output_directory,
      paste0(
        object_name,
        "_UMAP_diagnostic_summary"
      )
    ),
    width = 17,
    height = 14
  )
  
  p_elbow <- ElbowPlot(
    object,
    reduction = "pca",
    ndims = maximum_pcs
  ) +
    ggtitle(
      paste0(
        object_name,
        ": PCA elbow"
      )
    ) +
    theme_analysis()
  
  save_plot_both(
    p_elbow,
    file.path(
      output_directory,
      paste0(
        object_name,
        "_PCA_elbow"
      )
    ),
    width = 7,
    height = 5
  )
  
  logmsg(
    "%s: independent PCA and UMAP complete",
    object_name
  )
  
  object
}

# =============================================================================
# RUN INDEPENDENT ANALYSES
# =============================================================================

logmsg(
  "STEP 5: Running independent Yu2021 analysis"
)

yu <- process_independently(
  object = yu,
  object_name = "Yu2021",
  output_directory = yu_dir
)

logmsg(
  "STEP 6: Running independent internal analysis"
)

obj <- process_independently(
  object = obj,
  object_name = "Internal",
  output_directory = internal_dir
)

saveRDS(
  yu,
  file.path(
    yu_dir,
    "Yu2021_independent_PCA_UMAP.rds"
  )
)

saveRDS(
  obj,
  file.path(
    internal_dir,
    "Internal_independent_PCA_UMAP.rds"
  )
)

# =============================================================================
# PREPARE OBJECTS FOR COMBINED INTEGRATION
# =============================================================================

logmsg(
  "STEP 7: Preparing objects for combined integration"
)

# Retain all Yu healthy cells
yu_integration <- yu

# Exclude internal Spike cells from the primary integration
internal_cells_keep <- WhichCells(
  obj,
  expression =
    condition_common %in%
    CONFIG$internal_conditions_integrate
)

obj_integration <- subset(
  obj,
  cells = internal_cells_keep
)

cat("\nInternal integration conditions:\n")

print(
  table(
    obj_integration$condition_common,
    useNA = "ifany"
  )
)

# =============================================================================
# SHARED GENES
# =============================================================================

common_genes <- intersect(
  rownames(yu_integration),
  rownames(obj_integration)
)

if (length(common_genes) < 1000) {
  
  stop(
    "Too few common genes between datasets: ",
    length(common_genes)
  )
}

gene_overlap_summary <- tibble(
  metric = c(
    "Yu2021 genes",
    "Internal genes",
    "Common genes",
    "Yu2021-only genes",
    "Internal-only genes"
  ),
  value = c(
    nrow(yu_integration),
    nrow(obj_integration),
    length(common_genes),
    length(
      setdiff(
        rownames(yu_integration),
        rownames(obj_integration)
      )
    ),
    length(
      setdiff(
        rownames(obj_integration),
        rownames(yu_integration)
      )
    )
  )
)

print(
  gene_overlap_summary
)

write_csv(
  gene_overlap_summary,
  file.path(
    integration_dir,
    "gene_overlap_summary.csv"
  )
)

write_lines(
  common_genes,
  file.path(
    integration_dir,
    "common_genes.txt"
  )
)

yu_integration <- subset(
  yu_integration,
  features = common_genes
)

obj_integration <- subset(
  obj_integration,
  features = common_genes
)

# =============================================================================
# BALANCED INTEGRATION FEATURE SELECTION
#
# This selects features using both objects rather than allowing the much larger
# internal object to dominate feature selection.
# =============================================================================

logmsg(
  "STEP 8: Selecting balanced integration features"
)

integration_features <- SelectIntegrationFeatures(
  object.list = list(
    Yu2021 = yu_integration,
    Internal = obj_integration
  ),
  nfeatures =
    CONFIG$n_integration_features,
  assay = c(
    CONFIG$assay,
    CONFIG$assay
  )
)

integration_features <- intersect(
  integration_features,
  common_genes
)

write_lines(
  integration_features,
  file.path(
    integration_dir,
    "integration_features.txt"
  )
)

cat(
  "\nIntegration features:",
  length(integration_features),
  "\n"
)

# =============================================================================
# OPTIONAL ASSAY CLASS ALIGNMENT
#
# The Yu object may use a legacy Assay and the internal object may use Assay5.
# Convert Yu to Assay5 when supported.
# =============================================================================

if (
  inherits(
    obj_integration[[CONFIG$assay]],
    "Assay5"
  ) &&
  !inherits(
    yu_integration[[CONFIG$assay]],
    "Assay5"
  )
) {
  
  conversion_result <- try(
    {
      yu_integration[[CONFIG$assay]] <- as(
        yu_integration[[CONFIG$assay]],
        Class = "Assay5"
      )
    },
    silent = TRUE
  )
  
  if (
    inherits(
      conversion_result,
      "try-error"
    )
  ) {
    
    warning(
      "Yu RNA assay could not be converted to Assay5. ",
      "Attempting merge with the existing assay classes."
    )
  }
}

# =============================================================================
# MERGE OBJECTS
# =============================================================================

logmsg(
  "STEP 9: Merging Yu2021 and internal objects"
)

combined <- merge(
  x = yu_integration,
  y = obj_integration,
  merge.data = FALSE,
  project = "Yu2021_Internal_Duodenum"
)

DefaultAssay(combined) <- CONFIG$assay

# Join Seurat v5 layers after merge
if (
  inherits(
    combined[[CONFIG$assay]],
    "Assay5"
  ) &&
  "JoinLayers" %in%
  getNamespaceExports(
    "SeuratObject"
  )
) {
  
  combined <- SeuratObject::JoinLayers(
    combined,
    assay = CONFIG$assay
  )
}

# =============================================================================
# CHECK INTEGRATION DESIGN
# =============================================================================

integration_design <- combined@meta.data %>%
  count(
    dataset,
    condition_common,
    name = "n_cells"
  )

print(
  integration_design
)

write_csv(
  integration_design,
  file.path(
    integration_dir,
    "integration_design_cell_counts.csv"
  )
)

cat("\nDataset by condition:\n")

print(
  table(
    combined$dataset,
    combined$condition_common
  )
)

warning(
  paste0(
    "Dataset and biological condition are confounded: ",
    "Yu2021 contains Healthy cells, while Internal contains Treated/Untreated. ",
    "Use Harmony embeddings for visualization, neighborhood comparison, ",
    "or label alignment—not for differential expression."
  )
)

# =============================================================================
# COMBINED NORMALIZATION AND PCA
# =============================================================================

logmsg(
  "STEP 10: Running combined normalization and PCA"
)

combined <- NormalizeData(
  combined,
  assay = CONFIG$assay,
  normalization.method =
    CONFIG$normalization_method,
  scale.factor =
    CONFIG$scale_factor,
  verbose = FALSE
)

combined <- ScaleData(
  combined,
  assay = CONFIG$assay,
  features = integration_features,
  verbose = FALSE
)

maximum_combined_pcs <- min(
  CONFIG$n_pcs,
  length(integration_features) - 1,
  ncol(combined) - 1
)

combined <- RunPCA(
  combined,
  assay = CONFIG$assay,
  features = integration_features,
  npcs = maximum_combined_pcs,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 1234,
  verbose = FALSE
)

combined_dims <- CONFIG$dims_use[
  CONFIG$dims_use <= maximum_combined_pcs
]

# =============================================================================
# PRE-HARMONY UMAP
# =============================================================================

logmsg(
  "STEP 11: Running pre-Harmony UMAP"
)

combined <- RunUMAP(
  combined,
  reduction = "pca",
  dims = combined_dims,
  reduction.name = "umap.preharmony",
  reduction.key = "preUMAP_",
  seed.use = 1234,
  verbose = FALSE
)

p_pre_dataset <- DimPlot(
  combined,
  reduction = "umap.preharmony",
  group.by = "dataset",
  raster = CONFIG$raster
) +
  ggtitle(
    "Before Harmony: dataset"
  )

p_pre_condition <- DimPlot(
  combined,
  reduction = "umap.preharmony",
  group.by = "condition_common",
  raster = CONFIG$raster
) +
  ggtitle(
    "Before Harmony: condition"
  )

p_pre_celltype <- DimPlot(
  combined,
  reduction = "umap.preharmony",
  group.by = "celltype_L2_common",
  label = TRUE,
  repel = TRUE,
  raster = CONFIG$raster
) +
  ggtitle(
    "Before Harmony: PanGI L2"
  ) +
  NoLegend()

p_pre_technical <- DimPlot(
  combined,
  reduction = "umap.preharmony",
  group.by = "technical_sample",
  raster = CONFIG$raster
) +
  ggtitle(
    "Before Harmony: technical sample"
  )

p_pre_harmony <- (
  p_pre_dataset |
    p_pre_condition
) /
  (
    p_pre_celltype |
      p_pre_technical
  ) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    title = "Yu2021 and internal data before Harmony"
  ) &
  theme(
    legend.position = "bottom"
  )

save_plot_both(
  p_pre_harmony,
  file.path(
    integration_dir,
    "PreHarmony_UMAP_summary"
  ),
  width = 17,
  height = 13
)

# =============================================================================
# RUN HARMONY
#
# Only dataset is corrected:
#   Yu2021 versus Internal
#
# Not corrected:
#   technical_sample
#   biological_subject
#   donorID_unified
#   Subject_fixed
#   condition_common
#   Treatment
# =============================================================================

logmsg(
  "STEP 12: Running Harmony using dataset"
)

combined <- RunHarmony(
  object = combined,
  group.by.vars =
    CONFIG$harmony_group,
  reduction.use = "pca",
  dims.use = combined_dims,
  assay.use = CONFIG$assay,
  reduction.save = "harmony",
  project.dim = FALSE,
  plot_convergence = TRUE,
  verbose = TRUE
)

# =============================================================================
# HARMONY NEIGHBORS, CLUSTERS, AND UMAP
# =============================================================================

harmony_dimensions <- seq_len(
  min(
    max(combined_dims),
    ncol(
      Embeddings(
        combined,
        reduction = "harmony"
      )
    )
  )
)

combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = harmony_dimensions,
  graph.name = c(
    "harmony_nn",
    "harmony_snn"
  ),
  verbose = FALSE
)

combined <- FindClusters(
  combined,
  graph.name = "harmony_snn",
  resolution =
    CONFIG$clustering_resolution,
  cluster.name = "harmony_cluster",
  random.seed = 1234,
  verbose = FALSE
)

combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = harmony_dimensions,
  reduction.name = "umap.harmony",
  reduction.key = "hUMAP_",
  seed.use = 1234,
  verbose = FALSE
)

# =============================================================================
# POST-HARMONY DIAGNOSTIC PLOTS
# =============================================================================

logmsg(
  "STEP 13: Plotting post-Harmony results"
)

p_harmony_dataset <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "dataset",
  raster = CONFIG$raster
) +
  ggtitle(
    "After Harmony: dataset"
  )

p_harmony_condition <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "condition_common",
  raster = CONFIG$raster
) +
  ggtitle(
    "After Harmony: biological condition"
  )

p_harmony_celltype <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "celltype_L2_common",
  label = TRUE,
  repel = TRUE,
  raster = CONFIG$raster
) +
  ggtitle(
    "After Harmony: PanGI L2"
  ) +
  NoLegend()

p_harmony_technical <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "technical_sample",
  raster = CONFIG$raster
) +
  ggtitle(
    "After Harmony: technical sample"
  )

p_harmony_subject <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "biological_subject",
  raster = CONFIG$raster
) +
  ggtitle(
    "After Harmony: biological subject"
  )

p_harmony_cluster <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "harmony_cluster",
  label = TRUE,
  repel = TRUE,
  raster = CONFIG$raster
) +
  ggtitle(
    "After Harmony: graph clusters"
  ) +
  NoLegend()

p_harmony_summary <- (
  p_harmony_dataset |
    p_harmony_condition
) /
  (
    p_harmony_celltype |
      p_harmony_cluster
  ) /
  (
    p_harmony_technical |
      p_harmony_subject
  ) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    title = "Yu2021 and internal duodenum after dataset-level Harmony"
  ) &
  theme(
    legend.position = "bottom"
  )

save_plot_both(
  p_harmony_summary,
  file.path(
    integration_dir,
    "PostHarmony_UMAP_summary"
  ),
  width = 17,
  height = 18
)

# =============================================================================
# BEFORE VERSUS AFTER HARMONY
# =============================================================================

p_before_after_dataset <- (
  p_pre_dataset |
    p_harmony_dataset
) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    title = "Dataset mixing before and after Harmony"
  ) &
  theme(
    legend.position = "bottom"
  )

save_plot_both(
  p_before_after_dataset,
  file.path(
    integration_dir,
    "Dataset_before_after_Harmony"
  ),
  width = 16,
  height = 7
)

p_before_after_celltype <- (
  p_pre_celltype |
    p_harmony_celltype
) +
  plot_annotation(
    title = "PanGI L2 structure before and after Harmony"
  )

save_plot_both(
  p_before_after_celltype,
  file.path(
    integration_dir,
    "PanGI_L2_before_after_Harmony"
  ),
  width = 16,
  height = 7
)

# =============================================================================
# SPLIT HARMONY UMAP BY DATASET
# =============================================================================

p_harmony_split_dataset <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "celltype_L2_common",
  split.by = "dataset",
  raster = CONFIG$raster
) +
  plot_annotation(
    title = "Post-Harmony PanGI L2 structure split by dataset"
  )

save_plot_both(
  p_harmony_split_dataset,
  file.path(
    integration_dir,
    "PostHarmony_PanGI_L2_split_by_dataset"
  ),
  width = 18,
  height = 8
)

# =============================================================================
# CELL-TYPE MIXING SUMMARY
# =============================================================================

celltype_dataset_summary <- combined@meta.data %>%
  filter(
    !is.na(celltype_L2_common)
  ) %>%
  count(
    celltype_L2_common,
    dataset,
    name = "n_cells"
  ) %>%
  group_by(
    celltype_L2_common
  ) %>%
  mutate(
    fraction_within_celltype =
      n_cells /
      sum(n_cells)
  ) %>%
  ungroup()

write_csv(
  celltype_dataset_summary,
  file.path(
    integration_dir,
    "celltype_dataset_mixing_summary.csv"
  )
)

# =============================================================================
# SAVE FINAL OBJECTS
# =============================================================================

logmsg(
  "STEP 14: Saving processed objects"
)

saveRDS(
  yu,
  file.path(
    CONFIG$out_dir,
    "Yu2021_independent_processed.rds"
  )
)

saveRDS(
  obj,
  file.path(
    CONFIG$out_dir,
    "Internal_independent_processed.rds"
  )
)

saveRDS(
  combined,
  file.path(
    CONFIG$out_dir,
    "Yu2021_Internal_Harmony_dataset_corrected.rds"
  )
)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

final_summary <- bind_rows(
  
  tibble(
    object = "Yu2021 independent",
    n_cells = ncol(yu),
    n_genes = nrow(yu),
    n_technical_samples = n_distinct(
      yu$technical_sample,
      na.rm = TRUE
    ),
    n_biological_subjects = n_distinct(
      yu$biological_subject,
      na.rm = TRUE
    ),
    n_conditions = n_distinct(
      yu$condition_common,
      na.rm = TRUE
    )
  ),
  
  tibble(
    object = "Internal independent",
    n_cells = ncol(obj),
    n_genes = nrow(obj),
    n_technical_samples = n_distinct(
      obj$technical_sample,
      na.rm = TRUE
    ),
    n_biological_subjects = n_distinct(
      obj$biological_subject,
      na.rm = TRUE
    ),
    n_conditions = n_distinct(
      obj$condition_common,
      na.rm = TRUE
    )
  ),
  
  tibble(
    object = "Combined Harmony",
    n_cells = ncol(combined),
    n_genes = nrow(combined),
    n_technical_samples = n_distinct(
      combined$technical_sample,
      na.rm = TRUE
    ),
    n_biological_subjects = n_distinct(
      combined$biological_subject,
      na.rm = TRUE
    ),
    n_conditions = n_distinct(
      combined$condition_common,
      na.rm = TRUE
    )
  )
)

print(
  final_summary,
  width = Inf
)

write_csv(
  final_summary,
  file.path(
    CONFIG$out_dir,
    "final_analysis_summary.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    CONFIG$out_dir,
    "sessionInfo.txt"
  )
)

cat("\n")
cat("====================================================================\n")
cat("YU2021 + INTERNAL HARMONY PIPELINE COMPLETE\n")
cat("====================================================================\n")

cat("\nTechnical variables:\n")
cat("  Yu2021: sampleID\n")
cat("  Internal: orig.ident pooled GEX library\n")

cat("\nBiological subject variables:\n")
cat("  Yu2021: donorID_unified\n")
cat("  Internal: Subject_fixed\n")

cat("\nHarmony correction variable:\n")
cat("  dataset: Yu2021 versus Internal\n")

cat("\nNot corrected by Harmony:\n")
cat("  sampleID\n")
cat("  orig.ident\n")
cat("  donorID_unified\n")
cat("  Subject_fixed\n")
cat("  Treatment\n")

cat("\nInternal conditions included in integration:\n")
cat(
  "  ",
  paste(
    CONFIG$internal_conditions_integrate,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat("\nMain output object:\n")
cat(
  "  ",
  file.path(
    CONFIG$out_dir,
    "Yu2021_Internal_Harmony_dataset_corrected.rds"
  ),
  "\n",
  sep = ""
)

cat("\nImportant limitation:\n")
cat(
  "  Healthy status is specific to Yu2021, while treated/untreated status ",
  "is specific to the internal study.\n",
  sep = ""
)
cat(
  "  Use the Harmony embedding for visualization and alignment, not for ",
  "differential-expression testing.\n",
  sep = ""
)

cat("====================================================================\n")