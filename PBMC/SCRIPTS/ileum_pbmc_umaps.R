suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# =============================================================================
# INPUT
# =============================================================================

pbmc_file  <- "/gpfs/gsfs12/users/cheungf/ACEPD/PBMC/PBMC_RDS/object_final_annotated.rds"
ileum_file <- "/gpfs/gsfs12/users/cheungf/ACEPD/obj_lineage_Haemopoietic_cleaned_withDay.rds"

out_pdf <- "PBMC_ILEUM_UMAP_MonacoMain_similar.pdf"
out_png <- "PBMC_ILEUM_UMAP_MonacoMain_similar.png"

# =============================================================================
# LOAD
# =============================================================================

#pbmc <- readRDS(pbmc_file)
#ileum_cd45 <- readRDS(ileum_file)

# =============================================================================
# HELPERS
# =============================================================================

find_monaco_main_col <- function(obj) {
  candidates <- c(
    "Monaco_main",
    "monaco_main",
    "monaco_ann1",
    "monaco_main_annotation",
    "MonacoMain"
  )
  hit <- candidates[candidates %in% colnames(obj@meta.data)]
  if (length(hit) > 0) return(hit[1])
  
  stop(
    "Could not find Monaco main column. Available metadata columns:\n",
    paste(colnames(obj@meta.data), collapse = ", ")
  )
}

find_umap_name <- function(obj) {
  red_names <- Reductions(obj)
  
  if ("umap" %in% red_names) return("umap")
  
  hit <- red_names[str_detect(red_names, regex("umap", ignore_case = TRUE))]
  if (length(hit) > 0) return(hit[1])
  
  stop(
    "Could not find UMAP reduction. Available reductions:\n",
    paste(red_names, collapse = ", ")
  )
}

clean_labels <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("__", " ") %>%
    str_squish()
}

# =============================================================================
# DETECT COLUMNS
# =============================================================================

pbmc_monaco_col  <- find_monaco_main_col(pbmc)
ileum_monaco_col <- find_monaco_main_col(ileum_cd45)

pbmc_umap  <- find_umap_name(pbmc)
ileum_umap <- find_umap_name(ileum_cd45)

cat("PBMC Monaco column:", pbmc_monaco_col, "\n")
cat("Ileum Monaco column:", ileum_monaco_col, "\n")
cat("PBMC UMAP reduction:", pbmc_umap, "\n")
cat("Ileum UMAP reduction:", ileum_umap, "\n")

# =============================================================================
# BUILD PLOTTING DATA
# =============================================================================

pbmc_embed <- Embeddings(pbmc, reduction = pbmc_umap) %>%
  as.data.frame() %>%
  rownames_to_column("cell")

ileum_embed <- Embeddings(ileum_cd45, reduction = ileum_umap) %>%
  as.data.frame() %>%
  rownames_to_column("cell")

pbmc_df <- pbmc_embed %>%
  mutate(
    annotation = clean_labels(pbmc@meta.data[[pbmc_monaco_col]]),
    dataset = "PBMC"
  )

ileum_df <- ileum_embed %>%
  mutate(
    annotation = clean_labels(ileum_cd45@meta.data[[ileum_monaco_col]]),
    dataset = "Ileum CD45+"
  )

plot_df <- bind_rows(pbmc_df, ileum_df)

# standardize UMAP column names
umap_cols <- colnames(plot_df)[grepl("^UMAP_|^umap_", colnames(plot_df))]
if (length(umap_cols) < 2) {
  stop("Could not find at least 2 UMAP columns in embeddings.")
}

colnames(plot_df)[match(umap_cols[1], colnames(plot_df))] <- "UMAP_1"
colnames(plot_df)[match(umap_cols[2], colnames(plot_df))] <- "UMAP_2"

# keep same annotation levels across both panels
annotation_levels <- sort(unique(plot_df$annotation))
plot_df$annotation <- factor(plot_df$annotation, levels = annotation_levels)

# =============================================================================
# LABEL POSITIONS
# use median center per dataset x annotation
# =============================================================================

label_df <- plot_df %>%
  group_by(dataset, annotation) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n > 20)

# =============================================================================
# COLORS
# similar style to your figure
# =============================================================================

annotation_colors <- c(
  "B cells" = "#F8766D",
  "Basophils" = "#D89000",
  "CD4+ T cells" = "#A3A500",
  "CD8+ T cells" = "#39B600",
  "Dendritic cells" = "#00BF7D",
  "Monocytes" = "#00BFC4",
  "Neutrophils" = "#00B0F6",
  "NK cells" = "#9590FF",
  "Progenitors" = "#E76BF3",
  "T cells" = "#FF62BC"
)

# add fallback colors for any annotations not listed
missing_ann <- setdiff(annotation_levels, names(annotation_colors))
if (length(missing_ann) > 0) {
  extra_cols <- scales::hue_pal()(length(missing_ann))
  names(extra_cols) <- missing_ann
  annotation_colors <- c(annotation_colors, extra_cols)
}

annotation_colors <- annotation_colors[annotation_levels]

# =============================================================================
# PLOT
# =============================================================================

p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = annotation)) +
  geom_point(size = 0.8, alpha = 0.9) +
#  geom_text(
#    data = label_df,
#    aes(label = annotation),
#    color = "black",
#    size = 5,
#    fontface = "plain",
#    inherit.aes = FALSE
#  ) +
  facet_wrap(~ dataset, ncol = 2) +
  scale_color_manual(values = annotation_colors, drop = FALSE) +
  labs(
    title = "PBMC and Ileum CD45+ UMAPs colored by Monaco main annotation",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "annotation"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92")
  )

ggsave(out_pdf, p, width = 16, height = 11)
ggsave(out_png, p, width = 16, height = 11, dpi = 300)

p



combined <- merge(pbmc, ileum)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)

combined <- RunPCA(combined)
combined <- RunHarmony(combined, group.by.vars = "dataset")

combined <- RunUMAP(combined, reduction = "harmony")