# Robust UMAP IFNG feature plot with PanGI_L2 labels
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

# 0. prerequisites
DefaultAssay(obj) <- "RNA"   # set assay you want to use

# 1. UMAP embeddings (exists?)
umap_emb <- tryCatch(Embeddings(obj, "umap"), error = function(e) stop("UMAP embeddings not found. Run RunUMAP() first."))
umap_df <- as.data.frame(umap_emb)
# standardize column names to lowercase umap_1/umap_2 for downstream code
cn <- colnames(umap_df)
if (!any(tolower(cn) %in% c("umap_1","umap_2"))) {
  # try common variants
  if (all(c("UMAP_1","UMAP_2") %in% cn)) {
    umap_df <- umap_df %>% rename(umap_1 = UMAP_1, umap_2 = UMAP_2)
  } else if (all(c("UMAP_1","UMAP.1") %in% cn)) {
    umap_df <- umap_df %>% rename(umap_1 = UMAP_1, umap_2 = UMAP.2)
  } else {
    # fallback: take first two columns as UMAP coords (risky)
    warning("UMAP coordinate column names not standard; using first two columns as UMAP_1/UMAP_2")
    colnames(umap_df)[1:2] <- c("umap_1","umap_2")
  }
} else {
  # ensure lowercase names
  colnames(umap_df)[tolower(colnames(umap_df)) == "umap_1"] <- "umap_1"
  colnames(umap_df)[tolower(colnames(umap_df)) == "umap_2"] <- "umap_2"
}

umap_df$cell <- rownames(umap_df)

# 2. expression matrix: prefer slot 'data', fallback to 'counts' or 'scale.data'
assay_use <- DefaultAssay(obj)
expr_mat <- NULL
slot_used <- NULL
for (slot_name in c("data","counts","scale.data")) {
  m_try <- tryCatch(GetAssayData(obj, assay = assay_use, slot = slot_name), error = function(e) NULL)
  if (!is.null(m_try) && ncol(m_try) > 0 && nrow(m_try) > 0) {
    expr_mat <- m_try
    slot_used <- slot_name
    break
  }
}
if (is.null(expr_mat)) stop("No expression matrix found in slots 'data','counts','scale.data' for assay: ", assay_use)
message("Using assay '", assay_use, "' slot '", slot_used, "' with dims: ", paste(dim(expr_mat), collapse = " x "))

# 3. find IFNG row robustly
gene_candidates <- c("IFNG","Ifng","ifng")
gene_ifng <- intersect(gene_candidates, rownames(expr_mat))
if (length(gene_ifng) == 0) {
  # exact case-insensitive
  idx <- which(tolower(rownames(expr_mat)) == "ifng")
  if (length(idx) == 0) {
    # partial grep
    idx2 <- grep("ifng", tolower(rownames(expr_mat)), value = FALSE)
    if (length(idx2) == 0) {
      stop("IFNG not found in expression matrix rownames. Example rows: ", paste(head(rownames(expr_mat), 20), collapse = ", "))
    } else {
      gene_ifng <- rownames(expr_mat)[idx2[1]]
      message("Using IFNG-like match: ", gene_ifng)
    }
  } else {
    gene_ifng <- rownames(expr_mat)[idx[1]]
    message("Using IFNG match (case-insensitive): ", gene_ifng)
  }
} else {
  gene_ifng <- gene_ifng[1]
  message("Found IFNG as: ", gene_ifng)
}

# 4. align cell names (UMAP vs expression vs meta)
expr_cells <- colnames(expr_mat)
umap_cells  <- umap_df$cell
common_cells <- intersect(umap_cells, expr_cells)
if (length(common_cells) == 0) stop("No common cells between UMAP embeddings and expression matrix.")
message("Cells in common: ", length(common_cells), " / UMAP:", nrow(umap_df), " / expr:", ncol(expr_mat))

# subset and keep consistent order
umap_df <- umap_df[match(common_cells, umap_df$cell), , drop=FALSE]
expr_vec <- as.numeric(expr_mat[gene_ifng, common_cells])

# 5. attach expression and metadata
umap_df$IFNG <- expr_vec
# If PanGI_L2 exists in meta.data, attach it; else create NA
if ("PanGI_L2" %in% colnames(obj@meta.data)) {
  md_pan <- as.character(obj@meta.data[common_cells, "PanGI_L2"])
} else {
  md_pan <- rep(NA_character_, length(common_cells))
}
umap_df$PanGI_L2 <- md_pan
umap_df$Treatment <- if ("Treatment" %in% colnames(obj@meta.data)) as.character(obj@meta.data[common_cells, "Treatment"]) else NA_character_

# 6. centroids (median) and labels (filter small clusters)
centroids <- umap_df %>% group_by(PanGI_L2) %>% summarize(umap_1 = median(umap_1, na.rm=TRUE), umap_2 = median(umap_2, na.rm=TRUE), n = n()) %>% ungroup()
min_cells_for_label <- 20
centroids <- centroids %>% mutate(show = !is.na(PanGI_L2) & n >= min_cells_for_label)

# 7. plotting (order by IFNG so high expressers plotted last)
plot_df <- umap_df %>% arrange(IFNG)
p <- ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = IFNG), size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "D", na.value = "grey80", begin = 0.2, end = 0.95, name = paste0(gene_ifng, " (log)")) +
  theme_minimal(base_size = 12) +
  labs(title = paste0("UMAP: ", gene_ifng, " expression with PanGI_L2 labels")) +
  theme(legend.position = "right")

p <- p + geom_text(data = centroids %>% filter(show), aes(x = umap_1, y = umap_2, label = PanGI_L2), color = "black", size = 3.2, fontface = "bold") +
  geom_point(data = centroids %>% filter(show), aes(x = umap_1, y = umap_2), color = "white", shape = 21, fill = "black", size = 2, alpha = 0.6)

# save
ggsave("umap_IFNG_with_PanGI_L2_labels.png", p, width = 8, height = 6, dpi = 300)

# optional split by Treatment
if (!all(is.na(umap_df$Treatment))) {
  p_split <- ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = IFNG), size = 0.6, alpha = 0.8) +
    scale_color_viridis_c(option = "D", name = paste0(gene_ifng, " (log)")) +
    facet_wrap(~Treatment) + theme_minimal() + labs(title = paste0(gene_ifng, " by Treatment (UMAP)"))
  ggsave("umap_IFNG_by_Treatment.png", p_split, width = 10, height = 5, dpi = 300)
}

message("Done — saved umap_IFNG_with_PanGI_L2_labels.png (and split figure if Treatment present).")
