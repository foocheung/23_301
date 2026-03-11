library(Seurat)
library(ggplot2)
library(cowplot)

#-------------------------
# 1. Load objects
#-------------------------
obj_panGI <- readRDS("object_final_annotated_panGI_withDay.rds")
obj_heme  <- readRDS("obj_lineage_Haemopoietic_cleaned_withDay.rds")
obj_nonH  <- readRDS("obj_lineage_Non_Haemopoietic_cleaned_withDay.rds")

# Gene for CD45
gene <- "PTPRC"   # CD45

theme_set(theme_classic() + theme(legend.position = "none"))

#-------------------------
# 2. Function to plot CD45 with levels
#-------------------------
plot_CD45 <- function(obj, obj_name, reduction = "umap") {

  # check if gene exists
  if (!(gene %in% rownames(obj))) {
    stop(paste(gene, "not found in object:", obj_name))
  }

  # extract expression values
  vals <- FetchData(obj, vars = gene)[[1]]
  rng  <- range(vals, na.rm = TRUE)
  label_text <- paste0("CD45 (", gene, ") range: ",
                       round(rng[1], 2), " to ", round(rng[2], 2))

  FeaturePlot(
    obj,
    features = gene,
    reduction = reduction
  ) +
    ggtitle(paste0(obj_name, "\n", label_text)) +
    theme(plot.title = element_text(size = 14, face = "bold"))
}

#-------------------------
# 3. Generate the plots
#-------------------------
p1 <- plot_CD45(obj_panGI, "PanGI Annotated")
p2 <- plot_CD45(obj_heme, "Haemopoietic")
p3 <- plot_CD45(obj_nonH, "Non-Haemopoietic")

#-------------------------
# 4. Combine into one figure
#-------------------------
final_plot <- plot_grid(p1, p2, p3, ncol = 3)
final_plot

#-------------------------
# 5. Save PDF
#-------------------------
ggsave("CD45_FeaturePlots_panGI_heme_nonheme.pdf",
       final_plot,
       width = 18, height = 6)



##########


library(Seurat)
library(ggplot2)
library(cowplot)

#-------------------------
# 1. Load objects
#-------------------------
obj_panGI <- readRDS("object_final_annotated_panGI_withDay.rds")
obj_heme  <- readRDS("obj_lineage_Haemopoietic_cleaned_withDay.rds")
obj_nonH  <- readRDS("obj_lineage_Non_Haemopoietic_cleaned_withDay.rds")

# CD45 gene
gene_cd45 <- "PTPRC"   # CD45
meta_annotations <- c("HPCA_main", "Monaco_main", "PanGI_L1")
mito_col <- "percent.mito"    # assumes this column exists in meta.data



#-------------------------
# 2. Helper: CD45 FeaturePlot with range in title
#-------------------------
plot_CD45 <- function(obj, obj_name, reduction = "umap") {
  if (!(gene_cd45 %in% rownames(obj))) {
    warning(paste(gene_cd45, "not found in object:", obj_name))
    return(NULL)
  }
  vals <- FetchData(obj, vars = gene_cd45)[[1]]
  rng  <- range(vals, na.rm = TRUE)
  label_text <- paste0("CD45 (", gene_cd45, ") range: ",
                       round(rng[1], 2), " to ", round(rng[2], 2))
  
  FeaturePlot(
    obj,
    features = gene_cd45,
    reduction = reduction
  ) +
    ggtitle(paste0(obj_name, "\n", label_text)) +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

#-------------------------
# 3. Helper: percent.mito FeaturePlot
#-------------------------
plot_mito <- function(obj, obj_name, reduction = "umap") {
  if (!(mito_col %in% colnames(obj@meta.data))) {
    warning(paste(mito_col, "not found in meta.data for:", obj_name))
    return(NULL)
  }
  
  # Put percent.mito into a temporary feature slot for plotting
  obj[[mito_col]] <- obj@meta.data[[mito_col]]
  
  FeaturePlot(
    obj,
    features = mito_col,
    reduction = reduction
  ) +
    ggtitle(paste0(obj_name, " – ", mito_col)) +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

#-------------------------
# 4. Helper: UMAPs colored by metadata annotations
#-------------------------
plot_annotation_umaps <- function(obj, obj_name, reduction = "umap") {
  plots <- list()
  
  for (ann in meta_annotations) {
    if (ann %in% colnames(obj@meta.data)) {
      p <- DimPlot(
        obj,
        reduction = reduction,
        group.by = ann,
        label = TRUE,
        label.size = 3
      ) +
        ggtitle(paste0(obj_name, " – ", ann)) +
        theme(plot.title = element_text(size = 12, face = "bold"))
      plots[[ann]] <- p
    } else {
      warning(paste(ann, "not found in meta.data for:", obj_name))
    }
  }
  
  plots
}

#-------------------------
# 5. Wrapper: make all plots for one object
#-------------------------
make_all_plots_for_obj <- function(obj, obj_name) {
  ann_plots  <- plot_annotation_umaps(obj, obj_name)
  cd45_plot  <- plot_CD45(obj, obj_name)
  mito_plot  <- plot_mito(obj, obj_name)
  
  plot_list <- c(ann_plots,
                 list(CD45 = cd45_plot),
                 list(percent_mito = mito_plot))
  
  # Drop NULLs (missing columns/genes)
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
  
  # Arrange in a grid, 3 per row
  grid <- plot_grid(plotlist = plot_list, ncol = 3)
  grid
}

#-------------------------
# 6. Generate and save per-object PDFs
#-------------------------
p_panGI <- make_all_plots_for_obj(obj_panGI, "PanGI Annotated")
p_heme  <- make_all_plots_for_obj(obj_heme,  "Haemopoietic")
p_nonH  <- make_all_plots_for_obj(obj_nonH,  "Non-Haemopoietic")


#tapply(obj_nonH$CD45_expr, obj_nonH$PanGI_L2, median, na.rm = TRUE)
#table(obj_nonH$PanGI_L2)
#ggsave("UMAP_Annotations_CD45_percentMito_PanGI.pdf",
#       p_panGI, width = 18, height = 12)

immune_pangi <- c(
  "B_plasma",
  "Conventional_CD4",
  "Conventional_CD8",
  "Cycling_T/NK",
  "DC",
  "Macrophage",
  "Mature_B",
  "NK",
  "Treg",
  "Unconventional_T/ILC",
  "Lymphoid_stromal_cell"
)

# logical vector of cells to keep
keep_cells <- !(obj_nonH$PanGI_L2 %in% immune_pangi)

# subset by columns (cells)
obj_nonH <- obj_nonH[, keep_cells]



p_panGI <- make_all_plots_for_obj(obj_panGI, "PanGI Annotated")
p_heme  <- make_all_plots_for_obj(obj_heme,  "Haemopoietic")
p_nonH  <- make_all_plots_for_obj(obj_nonH,  "Non-Haemopoietic")


ggsave("UMAP_Annotations_CD45_percentMito_Haemopoietic.pdf",
       p_heme, width = 18, height = 12)

ggsave("UMAP_Annotations_CD45_percentMito_NonHaemopoietic.pdf",
       p_nonH, width = 18, height = 12)


###############3


# CD45 gene
gene_cd45 <- "PTPRC"   # CD45
#meta_annotations <- c("HPCA_main", "Monaco_main", "PanGI_L1")
meta_annotations <- c("Monaco_main", "PanGI_L1")
mito_col <- "percent.mito"    # assumes this column exists in meta.data

#-------------------------
# 2. Helper: CD45 FeaturePlot with range in title
#-------------------------
plot_CD45 <- function(obj, obj_name, reduction = "umap") {
  if (!(gene_cd45 %in% rownames(obj))) {
    warning(paste(gene_cd45, "not found in object:", obj_name))
    return(NULL)
  }
  vals <- FetchData(obj, vars = gene_cd45)[[1]]
  rng  <- range(vals, na.rm = TRUE)
  label_text <- paste0("CD45 (", gene_cd45, ") range: ",
                       round(rng[1], 2), " to ", round(rng[2], 2))
  
  FeaturePlot(
    obj,
    features = gene_cd45,
    reduction = reduction
  ) +
    ggtitle(paste0(obj_name, "\n", label_text)) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
}

#-------------------------
# 3. Helper: percent.mito FeaturePlot
#-------------------------
plot_mito <- function(obj, obj_name, reduction = "umap") {
  if (!(mito_col %in% colnames(obj@meta.data))) {
    warning(paste(mito_col, "not found in meta.data for:", obj_name))
    return(NULL)
  }
  
  # Put percent.mito into a temporary feature slot for plotting
  obj[[mito_col]] <- obj@meta.data[[mito_col]]
  
  FeaturePlot(
    obj,
    features = mito_col,
    reduction = reduction
  ) +
    ggtitle(paste0(obj_name, " – ", mito_col)) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
}

#-------------------------
# 4. Helper: UMAPs colored by metadata annotations
#-------------------------
plot_annotation_umaps <- function(obj, obj_name, reduction = "umap") {
  plots <- list()
  
  for (ann in meta_annotations) {
    if (ann %in% colnames(obj@meta.data)) {
      p <- DimPlot(
        obj,
        reduction = reduction,
        group.by = ann,
        label = TRUE,
        label.size = 3
      ) +
        ggtitle(paste0(obj_name, " – ", ann)) +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "none"
        )
      plots[[ann]] <- p
    } else {
      warning(paste(ann, "not found in meta.data for:", obj_name))
    }
  }
  
  plots
}

#-------------------------
# 5. Wrapper: make all plots for one object
#-------------------------
make_all_plots_for_obj <- function(obj, obj_name) {
  ann_plots  <- plot_annotation_umaps(obj, obj_name)
  cd45_plot  <- plot_CD45(obj, obj_name)
  mito_plot  <- plot_mito(obj, obj_name)
  
  plot_list <- c(ann_plots,
                 list(CD45 = cd45_plot),
                 list(percent_mito = mito_plot))
  
  # Drop NULLs (missing columns/genes)
  plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
  
  # Arrange in a grid, 3 per row
  grid <- plot_grid(plotlist = plot_list, ncol = 4)
  grid
}

#-------------------------
# 6. Generate per-object plots
#-------------------------
p_panGI <- make_all_plots_for_obj(obj_panGI, "PanGI Annotated")
p_heme  <- make_all_plots_for_obj(obj_heme,  "Haemopoietic")
p_nonH  <- make_all_plots_for_obj(obj_nonH,  "Non-Haemopoietic")

# (optional) Save to PDF if you want:
 ggsave("UMAP_Annotations_CD45_percentMito_PanGI.png", p_panGI, width = 16, height = 6)
 ggsave("UMAP_Annotations_CD45_percentMito_Haemopoietic.png", p_heme, width = 16, height = 6)
 ggsave("UMAP_Annotations_CD45_percentMito_NonHaemopoietic.png", p_nonH, width = 16, height = 6)
