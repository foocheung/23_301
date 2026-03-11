library(Seurat)
library(ggplot2)

## ------------------------------------
## Marker → cell-type mapping
## ------------------------------------
marker_info <- data.frame(
  gene = c(
    "DLL1",  "DEFA5", "PAX6", "MUC2",  "CHGA",   "SPDEF",  "NEUROD1",
    "ARX","DEFA6", "REG4",   "OLFM4", "LGR5",   "RBP2",   "FABP2",
    "ALPI",  "LYZ",   "BEST4",  "SPIB",  "TRPM5",  "TAS1R3"
  ),
  cell_type = c(
    "goblet",
    "Paneth",
    "enteroendocrine",
    "goblet",
    "enteroendocrine",
    "goblet / secretory",
    "enteroendocrine",
    "enteroendocrine",
    "Paneth",
    "deep crypt secretory",
    "stem / progenitor",
    "stem / progenitor",
    "enterocyte",
    "enterocyte",
    "enterocyte",
    "Paneth / antimicrobial",
    "BEST4+ enterocyte",
    "Paneth-like",
    "tuft / chemosensory",
    "tuft / chemosensory"
  ),
  stringsAsFactors = FALSE
)

marker_genes  <- marker_info$gene
marker_labels <- setNames(
  paste0(marker_info$gene, " – ", marker_info$cell_type),
  marker_info$gene
)

## ------------------------------------
## Violin plot: PanGI_L1 vs markers
## ------------------------------------
make_panGI_violin <- function(obj, obj_name,
                              group_col = "PanGI_L3",
                              assay = "RNA") {
  
  if (!(group_col %in% colnames(obj@meta.data))) {
    stop(paste("Column", group_col, "not found in meta.data for", obj_name))
  }
  
  # keep only markers that are present in this object
  feats <- intersect(marker_genes, rownames(obj))
  if (length(feats) == 0) {
    stop(paste("None of the marker genes found in", obj_name))
  }
  
  Idents(obj) <- obj[[group_col, drop = TRUE]]
  
  p_vln <- VlnPlot(
    obj,
    features = rev(feats),     # top-to-bottom order
    group.by = group_col,
    assay = assay,
    stack = TRUE,
    flip  = TRUE,
    pt.size = 0
  ) +
    ggtitle(paste0(obj_name, " – PanGI annotations vs epithelial markers")) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    ) +
    # relabel each gene strip as "GENE – cell type"
    scale_x_discrete(labels = marker_labels)
  
  p_vln
}

## Example: full PanGI object
p_panGI_violin <- make_panGI_violin(obj_nonH, "PanGI Annotated")
ggsave("3PanGI_PanGI_L1_violin_markers_annotated.png",
       p_panGI_violin, width = 12, height = 6, dpi = 300)


#######3

library(Seurat)
library(ggplot2)

## ------------------------------------
## Marker → cell-type mapping
## ------------------------------------
marker_info <- data.frame(
  gene = c(
    "DLL1",  "DEFA5", "PAX6", "MUC2",  "CHGA",   "SPDEF",  "NEUROD1",
    "ARX","DEFA6", "REG4",   "OLFM4", "LGR5",   "RBP2",   "FABP2",
    "ALPI",  "LYZ",   "BEST4",  "SPIB",  "TRPM5",  "TAS1R3"
  ),
  cell_type = c(
    "goblet",
    "Paneth",
    "enteroendocrine",
    "goblet",
    "enteroendocrine",
    "goblet / secretory",
    "enteroendocrine",
    "enteroendocrine",
    "Paneth",
    "deep crypt secretory",
    "stem / progenitor",
    "stem / progenitor",
    "enterocyte",
    "enterocyte",
    "enterocyte",
    "Paneth / antimicrobial",
    "BEST4+ enterocyte",
    "Paneth-like",
    "tuft / chemosensory",
    "tuft / chemosensory"
  ),
  stringsAsFactors = FALSE
)
marker_genes  <- marker_info$gene
marker_labels <- setNames(
  paste0(marker_info$gene, " – ", marker_info$cell_type),
  marker_info$gene
)

## ------------------------------------
## Violin plot: PanGI_L3 vs markers (with cell type labels)
## ------------------------------------
make_panGI_violin <- function(obj, obj_name,
                              group_col = "PanGI_L3",
                              assay = "RNA") {
  
  if (!(group_col %in% colnames(obj@meta.data))) {
    stop(paste("Column", group_col, "not found in meta.data for", obj_name))
  }
  
  # Keep only markers that are present in this object
  feats <- intersect(marker_genes, rownames(obj))
  if (length(feats) == 0) {
    stop(paste("None of the marker genes found in", obj_name))
  }
  
  Idents(obj) <- obj[[group_col, drop = TRUE]]
  
  # Create violin plot with NO individual cell points
  p_vln <- VlnPlot(
    obj,
    features = rev(feats),     # top-to-bottom order
    group.by = group_col,
    assay = assay,
    stack = TRUE,
    flip = TRUE,
    pt.size = 0                # NO points, just violins
  ) +
    ggtitle(paste0(obj_name, " – PanGI annotations vs epithelial markers")) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
    ) +
    # Y-axis labels show "GENE – cell type"
    scale_y_discrete(labels = marker_labels[rev(feats)])
  
  p_vln
}

## ------------------------------------
## Example usage
## ------------------------------------

# Generate violin plot with gene + cell type labels
p_panGI_violin <- make_panGI_violin(obj_nonH, "PanGI Annotated")

ggsave("3PanGI_PanGI_L3_violin_markers_annotated.png",
       p_panGI_violin, width = 12, height = 8, dpi = 300)



##############

library(Seurat)
library(ggplot2)
library(Seurat)
library(ggplot2)

## ------------------------------------
## Marker → cell-type mapping
## ------------------------------------
marker_info <- data.frame(
  gene = c(
    "DLL1",  "DEFA5", "PAX6", "MUC2",  "CHGA",   "SPDEF",  "NEUROD1",
    "ARX","DEFA6", "REG4",   "OLFM4", "LGR5",   "RBP2",   "FABP2",
    "ALPI",  "LYZ",   "BEST4",  "SPIB",  "TRPM5",  "TAS1R3"
  ),
  cell_type = c(
    "goblet",
    "Paneth",
    "enteroendocrine",
    "goblet",
    "enteroendocrine",
    "goblet / secretory",
    "enteroendocrine",
    "enteroendocrine",
    "Paneth",
    "deep crypt secretory",
    "stem / progenitor",
    "stem / progenitor",
    "enterocyte",
    "enterocyte",
    "enterocyte",
    "Paneth / antimicrobial",
    "BEST4+ enterocyte",
    "tuft/BEST4+",
    "tuft / chemosensory",
    "tuft / chemosensory"
  ),
  stringsAsFactors = FALSE
)

# hard-coded combined label column: "MARKER\ncell type"
marker_info$marker_label <- paste0(marker_info$gene, "\n", marker_info$cell_type)
marker_genes <- marker_info$gene

## ------------------------------------
## Violin plot: PanGI_L3 vs markers (marker + cell type on y-axis)
## ------------------------------------
make_panGI_violin <- function(obj,
                              obj_name = "PanGI Annotated",
                              group_col = "PanGI_L3",
                              assay = "RNA") {
  
  if (!(group_col %in% colnames(obj@meta.data))) {
    stop(paste("Column", group_col, "not found in meta.data for", obj_name))
  }
  
  DefaultAssay(obj) <- assay
  
  # keep only markers present in this object
  feats <- intersect(marker_genes, rownames(obj))
  if (length(feats) == 0) {
    stop(paste("None of the marker genes found in", obj_name))
  }
  
  # subset marker_info to match order of feats
  mi <- marker_info[match(feats, marker_info$gene), , drop = FALSE]
  
  # base Seurat violin plot (this gives you the “good” first plot)
  p <- VlnPlot(
    obj,
    features = feats,
    group.by = group_col,
    assay = assay,
    stack = TRUE,
    flip = TRUE,
    pt.size = 0
  ) +
    ggtitle(paste0(obj_name, " – PanGI_L3 vs epithelial markers")) +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1)
    )
  
  # --- KEY STEP: relabel the internal 'feature' column ---
  # p$data$feature currently holds just gene names; replace with "gene\ncelltype"
  p$data$feature <- factor(
    p$data$feature,
    levels = mi$gene,
    labels = mi$marker_label
  )
  
  # no extra facet_* call: VlnPlot already facets by `feature`
  # we just tweaked what `feature` looks like
  p
}

## ------------------------------------
## Example usage
## ------------------------------------
p_panGI_violin <- make_panGI_violin(obj_nonH, "Non-Haemopoietic")

ggsave(
  "PanGI_L3_violin_markers_marker_celltype.png",
  p_panGI_violin,
  width = 12,
  height = 12
 # dpi = 300
)
