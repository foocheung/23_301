library(Seurat)
library(future)

## 0) No futures, generous globals just in case
plan(sequential)
options(future.globals.maxSize = 32 * 1024^3)



qry  <- readRDS("../../MERGED_PIPELINE_OUTPUT_V5/object_final_annotated.rds")
ref <- readRDS("rds/healthy_reference/1_Healthy_Pan-GI_atlas_all_lineages_20241119_test_r2.rds")  # or lineage-specific file


# --- pick a label column that exists in the reference ---
colnames(ref@meta.data)[1:20]  # inspect


## 1) Make sure both objects are LogNormalized & have variable features
# (skip if you already did this)
ref <- NormalizeData(ref, verbose=FALSE)
qry <- NormalizeData(qry, verbose=FALSE)
if (length(VariableFeatures(ref)) == 0)  ref <- FindVariableFeatures(ref, nfeatures=3000, verbose=FALSE)
if (length(VariableFeatures(qry)) == 0)  qry <- FindVariableFeatures(qry, nfeatures=3000, verbose=FALSE)

## 2) Choose transfer features FIRST
features <- SelectIntegrationFeatures(object.list = list(ref, qry), nfeatures = 3000)

## 3) Trim to ONLY those features (this is the big win), drop heavy slots
ref2 <- DietSeurat(
  ref,
  assays     = "RNA",
  features   = features,   # <-- keep only transfer features
  counts     = FALSE,
  data       = TRUE,
  scale.data = FALSE,
  graphs     = NULL,
  dimreducs  = NULL        # we'll recompute PCA next
)
qry2 <- DietSeurat(
  qry,
  assays     = "RNA",
  features   = features,   # <-- match features exactly
  counts     = FALSE,
  data       = TRUE,
  scale.data = FALSE,
  graphs     = NULL,
  dimreducs  = NULL
)

## 4) Recompute PCA on trimmed objects (now tiny)
ref2 <- ScaleData(ref2, verbose=FALSE)
ref2 <- RunPCA(ref2, npcs=50, verbose=FALSE)
qry2 <- ScaleData(qry2, verbose=FALSE)
qry2 <- RunPCA(qry2, npcs=50, verbose=FALSE)

## 5) Anchors with fewer PCs
dims_use <- 1:30
anchors <- FindTransferAnchors(
  reference = ref2,
  query = qry2,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = dims_use,
  features = features,
  verbose = TRUE
)


##30 mins
#anchors <- FindTransferAnchors(
#  reference = ref2,
#  query = qry2,
#  normalization.method = "LogNormalize",
#  reduction = "rpca",                 # <-- faster path
# dims = dims_use,
#  features = features,
#  k.anchor = 3,       # default 5; smaller = faster (and still fine for mapping)
#  k.filter = 100,     # default 200; shrink neighbor filter
#  k.score  = 20,      # default 30; shrink scoring neighborhood
#  verbose = TRUE
#)


## 5) Transfer labels (example: L1/L2/L3)
pred_L1 <- TransferData(anchors, refdata = ref2[["level_1_annot", drop=TRUE]], dims = dims_use)
pred_L2 <- TransferData(anchors, refdata = ref2[["level_2_annot", drop=TRUE]], dims = dims_use)
pred_L3 <- TransferData(anchors, refdata = ref2[["level_3_annot", drop=TRUE]], dims = dims_use)

qry <- AddMetaData(qry, pred_L1)
qry <- AddMetaData(qry, pred_L2)
qry <- AddMetaData(qry, pred_L3)

cells <- colnames(qry)

# L1
pred_L1 <- as.data.frame(pred_L1)
rownames(pred_L1) <- if (is.null(rownames(pred_L1))) cells else rownames(pred_L1)
pred_L1 <- pred_L1[, 1:2, drop = FALSE]
colnames(pred_L1) <- c("PanGI_L1", "PanGI_L1_score")
qry <- AddMetaData(qry, pred_L1)

# L2
pred_L2 <- as.data.frame(pred_L2)
rownames(pred_L2) <- if (is.null(rownames(pred_L2))) cells else rownames(pred_L2)
pred_L2 <- pred_L2[, 1:2, drop = FALSE]
colnames(pred_L2) <- c("PanGI_L2", "PanGI_L2_score")
qry <- AddMetaData(qry, pred_L2)

# L3
pred_L3 <- as.data.frame(pred_L3)
rownames(pred_L3) <- if (is.null(rownames(pred_L3))) cells else rownames(pred_L3)
pred_L3 <- pred_L3[, 1:2, drop = FALSE]
colnames(pred_L3) <- c("PanGI_L3", "PanGI_L3_score")
qry <- AddMetaData(qry, pred_L3)



p1 <- DimPlot(
  qry, group.by = "PanGI_L2", reduction = "umap", label = TRUE
) +
  theme_void()+
  ggtitle("Pan-GI L2 annotations")

# CD45 gene name = PTPRC in the Pan-GI atlas
p2 <- FeaturePlot(
  qry, features = "PTPRC", reduction = "umap"
) +
  scale_color_viridis_c() +
  ggtitle("CD45 (PTPRC) expression")

# Combine plots
p1 + p2



#############


#obj  <- readRDS("../../obj_singletons_filtered.rds")
#obj  <- readRDS("../../MERGED_PIPELINE_OUTPUT_V4/object_final_annotated.rds")
obj  <- readRDS("../../MERGED_PIPELINE_OUTPUT_V5/object_final_annotated.rds")
# L1
pred_L1 <- as.data.frame(pred_L1)
rownames(pred_L1) <- if (is.null(rownames(pred_L1))) cells else rownames(pred_L1)
pred_L1 <- pred_L1[, 1:2, drop = FALSE]
colnames(pred_L1) <- c("PanGI_L1", "PanGI_L1_score")
obj <- AddMetaData(obj, pred_L1)

# L2
pred_L2 <- as.data.frame(pred_L2)
rownames(pred_L2) <- if (is.null(rownames(pred_L2))) cells else rownames(pred_L2)
pred_L2 <- pred_L2[, 1:2, drop = FALSE]
colnames(pred_L2) <- c("PanGI_L2", "PanGI_L2_score")
obj <- AddMetaData(obj, pred_L2)

# L3
pred_L3 <- as.data.frame(pred_L3)
rownames(pred_L3) <- if (is.null(rownames(pred_L3))) cells else rownames(pred_L3)
pred_L3 <- pred_L3[, 1:2, drop = FALSE]
colnames(pred_L3) <- c("PanGI_L3", "PanGI_L3_score")
obj <- AddMetaData(obj, pred_L3)


obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
obj@meta.data %>% colnames()

saveRDS(obj, "../../MERGED_PIPELINE_OUTPUT_V5/object_final_annotated_panGI.rds")


library(ggplot2)
library(patchwork)

p_pan <- DimPlot(obj, reduction = "umap", group.by = "PanGI_L2",
                 label = TRUE, repel = TRUE) + ggtitle("Pan-GI (L2)") + NoLegend()

p_tlas <- DimPlot(obj, reduction = "umap", group.by = "HPCA_main",
                 label = TRUE, repel = TRUE) + ggtitle("HPCA_main") + NoLegend()

p_mon <- DimPlot(obj, reduction = "umap", group.by = "Monaco_main",
                  label = TRUE, repel = TRUE) + ggtitle("Monaco_main") + NoLegend()

p_pan + p_mon + p_tlas


p_pan <- DimPlot(obj, reduction = "umap", group.by = "PanGI_L1",
                 label = TRUE, repel = TRUE) + ggtitle("Pan-GI (L1)") + NoLegend()

p_pan2 <- DimPlot(obj, reduction = "umap", group.by = "PanGI_L2",
                 label = TRUE, repel = TRUE) + ggtitle("Pan-GI (L2)") + NoLegend()


p_pan3 <- DimPlot(obj, reduction = "umap", group.by = "PanGI_L3",
                 label = TRUE, repel = TRUE) + ggtitle("Pan-GI (L3)") + NoLegend()


p_pan + p_pan2 + p_pan3



##################


# 1) Choose gene & assay/slot (CD45 = PTPRC)
gene  <- "PTPRC"
assay <- "RNA"      # change if your expression is in another assay
slot  <- "data"     # "data" = log-normalized; use "counts" if you want raw

stopifnot(gene %in% rownames(GetAssayData(obj, assay = assay)))

# 2) Pull expression vector and combine with L1 labels
expr_vec <- GetAssayData(obj, assay = assay, slot = slot)[gene, , drop = TRUE]
df <- data.frame(
  PanGI_L1 = obj$PanGI_L1,
  CD45 = as.numeric(expr_vec),
  check.names = FALSE
)

# 3) Median per L1 (plus N for context)
med_cd45 <- df %>%
  group_by(PanGI_L1) %>%
  summarise(
    median_CD45 = median(CD45, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_CD45))



med_cd45
# A tibble: 7 × 3
PanGI_L1       median_CD45     n
<chr>                <dbl> <int>
1 T and NK cells        3.15 59060
2 Myeloid               2.14  2789
3 B and B plasma        2.02 13172
4 Endothelial           0      253
5 Epithelial            0    64494
6 Mesenchymal           0      541
7 Neural                0       15




qry2 <- MapQuery(anchors, qry2, ref2, refdata=NULL,
                 reference.reduction="pca", reduction.model="umap",
                 new.reduction.name="ref.umap")
