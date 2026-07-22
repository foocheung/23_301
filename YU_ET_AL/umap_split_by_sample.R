# =============================================================================
# POST-HARMONY UMAP:
# PAN-GI L2 ANNOTATION SPLIT BY TECHNICAL SAMPLE
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# Load the integrated Harmony object if it is not already in memory
#combined <- readRDS(
#  "YU2021_INTERNAL_HARMONY/Yu2021_Internal_Harmony_dataset_corrected.rds"
#)

# Check that the needed reduction and metadata exist
stopifnot(
  "umap.harmony" %in% Reductions(combined),
  "celltype_L2_common" %in% colnames(combined@meta.data),
  "technical_sample" %in% colnames(combined@meta.data)
)

# =============================================================================
# 1. ANNOTATION SPLIT BY TECHNICAL SAMPLE
# =============================================================================

p_annotation_by_sample <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "celltype_L2_common",
  split.by = "technical_sample",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  ncol = 4
) +
  plot_annotation(
    title = "Post-Harmony PanGI L2 annotation split by technical sample"
  )

print(p_annotation_by_sample)

ggsave(
  filename = "YU2021_INTERNAL_HARMONY/04_HARMONY_INTEGRATION/PostHarmony_PanGI_L2_split_by_technical_sample.png",
  plot = p_annotation_by_sample,
  width = 22,
  height = 18,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "YU2021_INTERNAL_HARMONY/04_HARMONY_INTEGRATION/PostHarmony_PanGI_L2_split_by_technical_sample.pdf",
  plot = p_annotation_by_sample,
  width = 22,
  height = 18,
  bg = "white"
)

# =============================================================================
# 2. ANNOTATION SPLIT BY DATASET
# =============================================================================

p_annotation_by_dataset <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "celltype_L2_common",
  split.by = "dataset",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  ncol = 2
) +
  plot_annotation(
    title = "Post-Harmony PanGI L2 annotation split by dataset"
  )

print(p_annotation_by_dataset)

ggsave(
  filename = "YU2021_INTERNAL_HARMONY/04_HARMONY_INTEGRATION/PostHarmony_PanGI_L2_split_by_dataset.png",
  plot = p_annotation_by_dataset,
  width = 18,
  height = 8,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# 3. ANNOTATION SPLIT BY BIOLOGICAL SUBJECT
# =============================================================================

p_annotation_by_subject <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "celltype_L2_common",
  split.by = "biological_subject",
  label = FALSE,
  raster = TRUE,
  ncol = 4
) +
  plot_annotation(
    title = "Post-Harmony PanGI L2 annotation split by biological subject"
  )

print(p_annotation_by_subject)

ggsave(
  filename = "YU2021_INTERNAL_HARMONY/04_HARMONY_INTEGRATION/PostHarmony_PanGI_L2_split_by_subject.png",
  plot = p_annotation_by_subject,
  width = 22,
  height = 18,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# 4. SAMPLE IDENTITY ON THE POST-HARMONY UMAP
# =============================================================================

p_sample_identity <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "technical_sample",
  raster = TRUE
) +
  ggtitle("Post-Harmony UMAP by technical sample")

print(p_sample_identity)

ggsave(
  filename = "YU2021_INTERNAL_HARMONY/04_HARMONY_INTEGRATION/PostHarmony_UMAP_by_technical_sample.png",
  plot = p_sample_identity,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# 5. PAN-GI L2 ANNOTATION WITH SAMPLE PANELS, NO LABELS
#
# This is often easier to read when there are many sample panels.
# =============================================================================

p_annotation_by_sample_clean <- DimPlot(
  combined,
  reduction = "umap.harmony",
  group.by = "celltype_L2_common",
  split.by = "technical_sample",
  label = FALSE,
  raster = TRUE,
  ncol = 4
) +
  plot_annotation(
    title = "Post-Harmony PanGI L2 annotation split by technical sample"
  ) &
  theme(
    legend.position = "bottom"
  )

#print(p_annotation_by_sample_clean)

ggsave(
  filename = "YU2021_INTERNAL_HARMONY/04_HARMONY_INTEGRATION/PostHarmony_PanGI_L2_split_by_technical_sample_clean.png",
  plot = p_annotation_by_sample_clean,
  width = 22,
  height = 18,
  dpi = 300,
  bg = "white"
)