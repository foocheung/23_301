suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# =============================================================================
# LOAD OBJECT
# =============================================================================
COMBINED_RDS <- "PBMC_RDS/object_final_annotated.rds"

obj <- readRDS(COMBINED_RDS)

# =============================================================================
# CHECK METADATA COLUMNS (IMPORTANT)
# =============================================================================
colnames(obj@meta.data)

# You likely have:
# seurat_clusters
# Monaco_fine_pruned (or similar)

# Rename here if needed
CLUSTER_COL <- "cluster"
ANNOT_COL   <- "Monaco_fine_pruned"

# =============================================================================
# CLEAN ANNOTATION NAMES
# =============================================================================
clean_names <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_]", "")
}

obj@meta.data <- obj@meta.data %>%
  mutate(
    cluster = as.character(.data[[CLUSTER_COL]]),
    celltype = clean_names(.data[[ANNOT_COL]])
  )

# =============================================================================
# 1. GET CELLTYPES_KEEP VECTOR
# =============================================================================
CELLTYPES_KEEP <- obj@meta.data %>%
  pull(celltype) %>%
  unique() %>%
  sort()

# print in copy-paste format
cat(
  "CELLTYPES_KEEP <- c(\n",
  paste0('  "', CELLTYPES_KEEP, '"', collapse = ",\n"),
  "\n)\n"
)

# =============================================================================
# 2. CLUSTER ↔ ANNOTATION TABLE
# =============================================================================
cluster_annotation_table <- obj@meta.data %>%
  count(cluster, celltype) %>%
  group_by(cluster) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(cluster, desc(freq))

# save
write.csv(cluster_annotation_table,
          "cluster_vs_monaco_fine_table.csv",
          row.names = FALSE)

# =============================================================================
# 3. DOMINANT CELLTYPE PER CLUSTER
# =============================================================================
cluster_dominant <- cluster_annotation_table %>%
  group_by(cluster) %>%
  slice_max(freq, n = 1) %>%
  ungroup()

write.csv(cluster_dominant,
          "cluster_dominant_celltype.csv",
          row.names = FALSE)

# =============================================================================
# 4. OPTIONAL: VISUALIZE MIXING (VERY USEFUL)
# =============================================================================
p <- ggplot(cluster_annotation_table,
            aes(x = cluster, y = freq, fill = celltype)) +
  geom_bar(stat = "identity") +
  
  # 🔥 LABEL WITH CELLTYPE NAMES
  geom_text(
    aes(label = ifelse(freq > 0.08, celltype, "")),  # only show big enough segments
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white",angle = 90     
  ) +
  
  theme_bw() +
  labs(
    title = "Cluster → Monaco Fine Cell Type Composition",
    y = "Fraction",
    x = "Seurat Cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("cluster_annotation_labeled_celltypes.png", p, width = 12, height = 7)


ggsave("cluster_annotation_mixing_plot.png", p, width = 10, height = 6)

# =============================================================================
# 5. OPTIONAL: QUICK UMAP CHECK
# =============================================================================
p1 <- DimPlot(obj, group.by = CLUSTER_COL, label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(obj, group.by = ANNOT_COL, label = TRUE) + ggtitle("Monaco Fine")

ggsave("umap_cluster_vs_annotation.png", p1 + p2, width = 12, height = 6)

# =============================================================================
# DONE
# =============================================================================
cat("Done.\n")



p <- ggplot(cluster_annotation_table,
            aes(x = cluster, y = freq, fill = celltype)) +
  geom_bar(stat = "identity") +

  # 🔥 ADD LABELS INSIDE
  geom_text(
    aes(label = ifelse(freq > 0.05, round(freq, 2), "")),  # only label bigger segments
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  ) +

  theme_bw() +
  labs(
    title = "Cluster → Cell Type Composition",
    y = "Fraction",
    x = "Seurat Cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("cluster_annotation_mixing_plot_labeled.png", p, width = 10, height = 6)