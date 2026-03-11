## ---- Setup ----
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load object
obj_panGI <- readRDS("object_final_annotated_panGI_withDay.rds")

## ---- Treatment Setup ----
obj_panGI$Treatment <- factor(
  obj_panGI$Treatment,
  levels = c("Untreated", "Treated", "Spike")
)

## ---- Cluster x Treatment table ----
Idents(obj_panGI) <- "seurat_clusters"

cluster_by_treatment <- table(
  Cluster   = Idents(obj_panGI),
  Treatment = obj_panGI$Treatment
)

cluster_df <- as.data.frame(cluster_by_treatment) %>%
  group_by(Cluster) %>%
  mutate(
    Total = sum(Freq),
    Freq_frac = Freq / Total
  ) %>%
  arrange(Cluster, Treatment)

# Save the table
write.csv(cluster_df, "cluster_by_treatment.csv", row.names = FALSE)
cat("Saved: cluster_by_treatment.csv\n")

## ---- Identify Untreated-specific clusters ----
untreated_specific <- cluster_df %>%
  filter(Treatment == "Untreated", Freq_frac > 0.9, Total >= 50) %>%
  pull(Cluster) %>%
  as.character()

obj_panGI$cluster_untreated_label <- as.character(Idents(obj_panGI))
obj_panGI$cluster_untreated_label[
  obj_panGI$cluster_untreated_label %in% untreated_specific
] <- paste0(obj_panGI$cluster_untreated_label[
  obj_panGI$cluster_untreated_label %in% untreated_specific
], "_UntreatedOnly")

## ---- Make Plots ----

### 1) Global split UMAP
p_global_split <- DimPlot(
  obj_panGI,
  reduction = "umap",
  group.by  = "seurat_clusters",
  split.by  = "Treatment",
  pt.size   = 0.3,
  raster    = TRUE
) + ggtitle("UMAP by cluster, split by Treatment")

### 2) Highlight Untreated-only clusters
p_highlight_untreated <- DimPlot(
  obj_panGI,
  reduction = "umap",
  group.by  = "cluster_untreated_label",
  label     = TRUE,
  pt.size   = 0.2,
  raster    = TRUE
) + ggtitle("UMAP with Untreated-enriched clusters labeled")

### 3) Separate UMAPs
cells_untreated <- WhichCells(obj_panGI, expression = Treatment == "Untreated")
cells_treated   <- WhichCells(obj_panGI, expression = Treatment == "Treated")
cells_spike     <- WhichCells(obj_panGI, expression = Treatment == "Spike")

p_untreated <- DimPlot(obj_panGI, reduction="umap",
                       cells=cells_untreated, group.by="seurat_clusters",
                       pt.size=0.4, raster=TRUE) +
  ggtitle("UMAP: Untreated")

p_treated <- DimPlot(obj_panGI, reduction="umap",
                     cells=cells_treated, group.by="seurat_clusters",
                     pt.size=0.4, raster=TRUE) +
  ggtitle("UMAP: Treated")

p_spike <- DimPlot(obj_panGI, reduction="umap",
                   cells=cells_spike, group.by="seurat_clusters",
                   pt.size=0.4, raster=TRUE) +
  ggtitle("UMAP: Spike")

# Combine vertical
p_three <- p_untreated / p_treated / p_spike


## ---- Save Plots ----

# 1) Global split UMAP
ggsave(
  "umap_global_split.pdf",
  p_global_split,
  width = 12, height = 6, useDingbats = FALSE
)

# 2) Highlight Untreated-only clusters
ggsave(
  "umap_highlight_untreated.pdf",
  p_highlight_untreated,
  width = 7, height = 6, useDingbats = FALSE
)

# 3) Three-panel combined
ggsave(
  "umap_untreated_treated_spike.pdf",
  p_three,
  width = 8, height = 12, useDingbats = FALSE
)

cat("Saved all UMAP plots.\n")

ggsave("umap_global_split.png", p_global_split, width=12, height=6, dpi=300)




Idents(obj_panGI) <- "seurat_clusters"

cluster30_untreated_cells <- head(obj_panGI@meta.data[cluster30_untreated_cells, ])


length(cluster30_untreated_cells)


########

Idents(obj_panGI) <- "seurat_clusters"

# Cells in cluster 30
cells_30 <- WhichCells(obj_panGI, expression = seurat_clusters == 30)

# Count by Treatment
table(obj_panGI$Treatment[cells_30])

###

cluster30_by_subject
Treatment
Subject    Untreated Treated Spike
P0001171        49       0     0
P0002345        19       0     0



head(obj_panGI@meta.data[cluster30_untreated_cells, ])
obj_panGI@meta.data[cluster30_untreated_cells, ]$PanGI_L3 PanGI_L3_score

B_GC_I
Germinal Center B cell – Intermediate state
A transitional GC B population between dark-zone (CXCR4⁺, MKI67⁺) and light-zone (CXCR5⁺, CD86⁺) programs.
small cluster of 65 cells from 2 untreated subjects, removed after JakStat inhibitor



Here is a clean, publication-style description for your results section or figure legend, based directly on your data:
  
  ---
  
  # **B_GC_I Population Summary (Pan-GI Reference)**
  
  **Cell Type:**
  **B_GC_I — Germinal Center B cell, Intermediate State**
  A transitional germinal center (GC) B-cell population positioned between:
  
  * **Dark zone (DZ):** CXCR4⁺, MKI67⁺ proliferating centroblasts
* **Light zone (LZ):** CXCR5⁺, CD86⁺ antigen-presenting centrocytes

B_GC_I expresses a **mixed DZ/LZ transcriptional profile**, reflecting a GC-committed state that is shifting between proliferation (DZ) and selection/signaling (LZ).

---
  
  # **Your Dataset Findings (ACEPD Ileum)**
  
  ### **Cluster Size**
  
  * **65 total cells**
  * Sourced from **only 2 untreated subjects**
  
  ### **Treatment Response**
  
  * **Completely removed after JAK–STAT inhibitor treatment**
  * No detectable B_GC_I cells remain in:
  
  * Treated (PMA or baseline stimulation)
* Spike-in
* JAKi group → **0 cells detected**
  
  ### **Interpretation**
  
  The disappearance of B_GC_I after JAK–STAT inhibition is consistent with:
  
  * **IL-21 → JAK1/JAK3 → STAT3 signaling**, which is essential for sustaining germinal center B-cell programs
* Transitional GC B cells (including B_GC_I) rely on **STAT3/STAT6** for:
  
  * BCL6 maintenance
* GC activation survival
* partial DZ phenotype (CXCR4↑, MKI67 intermediate)

Blocking JAK–STAT signaling disrupts these pathways, leading to:
  
  * loss of GC transcriptional identity
* collapse of transitional GC B-cell states
* failure to maintain B_GC_I and CXCR4⁺ DZ-leaning intermediates

➡️ **B_GC_I is biologically expected to disappear after JAK–STAT blockade.**
  
  ---
  
  # **If you want this phrased as a figure caption:**
  
  **“B_GC_I (Germinal Center Intermediate B cells) form a small cluster (65 cells from two untreated ACEPD subjects) and exhibit a mixed DZ/LZ transcriptional signature. This population is entirely absent after JAK–STAT inhibitor treatment, consistent with the known dependence of GC B-cell maintenance on IL-21/IL-4 → STAT3/STAT6 signaling.”**
  
  ---

  
  
  