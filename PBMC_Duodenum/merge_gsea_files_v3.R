##CORRECT FO DUODENUM
#for the cluster name make shure to sub space with an underscore eg "_"
#Vd2 gd T cells
#Vd2_gd_T_cells

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(fs)
})

# =============================================================================
# PATHS
# =============================================================================
###STILL NEEED TO MANUALLY FIX CELL NAMES TO MAKE IT CONSISTENT!!!!!

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(fs)
})

# =============================================================================
# PATHS
# =============================================================================
PBMC_ROOT <- "deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned/paired_cross_pool/gsea_by_celltype"
ILEUM_ROOT <- "/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/OUTPUT_20251128_141815_WITH_BATCH/Monaco_fine/gsea_results"

OUT_FILE <- "PBMC_ILEUM_COMBINED_GSEA_ALL_FINE_ANN.csv"

# =============================================================================
# HELPERS
# =============================================================================
clean_cluster_name <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_")
}

extract_pbmc_celltype <- function(f) {
  # safest: use parent folder name under gsea_by_celltype
  parts <- str_split(path_norm(f), "/")[[1]]
  idx <- which(parts == "gsea_by_celltype")
  
  if (length(idx) > 0 && idx[length(idx)] < length(parts)) {
    return(clean_cluster_name(parts[idx[length(idx)] + 1]))
  }
  
  # fallback to filename
  fname <- basename(f)
  out <- fname %>%
    str_remove("_MildTreated_vs_MildUntreated_paired_.*$") %>%
    clean_cluster_name()
  
  out
}

extract_pbmc_geneset <- function(fname) {
  case_when(
    str_detect(fname, "GO_BP") ~ "GO:BP",
    str_detect(fname, "Reactome") ~ "REACTOME",
    str_detect(fname, "Hallmark") ~ "HALLMARK",
    TRUE ~ "UNKNOWN"
  )
}

extract_ileum_cluster <- function(fname) {
  out <- fname %>%
    str_extract("cluster_.*?_GSEA") %>%
    str_remove("^cluster_") %>%
    str_remove("_GSEA$") %>%
    clean_cluster_name()
  
  out
}

extract_ileum_geneset <- function(fname) {
  case_when(
    str_detect(fname, "GO") ~ "GO:BP",
    str_detect(fname, "REACTOME") ~ "REACTOME",
    str_detect(fname, "HALLMARK") ~ "HALLMARK",
    TRUE ~ "UNKNOWN"
  )
}

# =============================================================================
# 1. LOAD PBMC FILES (UNFILTERED ONLY)
# =============================================================================
pbmc_files <- dir_ls(PBMC_ROOT, recurse = TRUE, regexp = "_results\\.csv$")
pbmc_files <- pbmc_files[!str_detect(pbmc_files, "_sig\\.csv")]

load_pbmc <- function(f) {
  df <- fread(f)
  fname <- basename(f)
  
  celltype <- extract_pbmc_celltype(f)
  geneset  <- extract_pbmc_geneset(fname)
  
  df$cluster <- celltype
  df$geneset <- geneset
  df$annotation <- "Monaco_fine"
  df$dataset <- "PBMC"
  
  # add missing ileum-style columns if absent
  if (!"leading_edge_symbols" %in% colnames(df)) df$leading_edge_symbols <- NA_character_
  if (!"leading_edge_entrez" %in% colnames(df)) df$leading_edge_entrez <- NA_character_
  
  df
}

pbmc_list <- map(pbmc_files, load_pbmc)
pbmc_all <- bind_rows(pbmc_list)

# =============================================================================
# 2. LOAD ILEUM FILES
# =============================================================================
ileum_files <- dir_ls(ILEUM_ROOT, recurse = TRUE, regexp = "cluster_.*Treated_vs_Untreated\\.csv$")

load_ileum <- function(f) {
  df <- fread(f)
  fname <- basename(f)
  
  cluster <- extract_ileum_cluster(fname)
  geneset <- extract_ileum_geneset(fname)
  
  df$cluster <- cluster
  df$geneset <- geneset
  df$annotation <- "Monaco_fine"
  df$dataset <- "ILEUM"
  
  # add missing PBMC-style columns if absent
  if (!"leading_edge_symbols" %in% colnames(df)) df$leading_edge_symbols <- NA_character_
  if (!"leading_edge_entrez" %in% colnames(df)) df$leading_edge_entrez <- NA_character_
  if (!"leading_edge" %in% colnames(df)) df$leading_edge <- NA_character_
  
  df
}

ileum_list <- map(ileum_files, load_ileum)
ileum_all <- bind_rows(ileum_list)

# =============================================================================
# 3. STANDARDIZE COLUMN ORDER
# =============================================================================
common_cols <- c(
  "ID", "Description", "setSize", "enrichmentScore", "NES",
  "pvalue", "p.adjust", "qvalue", "rank",
  "leading_edge", "leading_edge_symbols", "leading_edge_entrez",
  "cluster", "geneset", "annotation", "dataset"
)

# add any missing columns before subsetting
for (col in common_cols) {
  if (!col %in% colnames(pbmc_all)) pbmc_all[[col]] <- NA
  if (!col %in% colnames(ileum_all)) ileum_all[[col]] <- NA
}

pbmc_all <- pbmc_all[, common_cols, with = FALSE]
ileum_all <- ileum_all[, common_cols, with = FALSE]

# =============================================================================
# 4. MERGE
# =============================================================================
combined <- bind_rows(pbmc_all, ileum_all)

# =============================================================================
# 5. ADD UNIQUE IDENTIFIERS
# =============================================================================
combined <- combined %>%
  mutate(
    cluster = clean_cluster_name(cluster),
    pathway_id = paste(ID, Description, sep = "|"),
    celltype_dataset = paste(dataset, cluster, sep = "_")
  )

# =============================================================================
# 6. SAVE
# =============================================================================
fwrite(combined, OUT_FILE)

cat("Saved:", OUT_FILE, "\n")
cat("PBMC rows:", nrow(pbmc_all), "\n")
cat("ILEUM rows:", nrow(ileum_all), "\n")
cat("Combined rows:", nrow(combined), "\n")










#########################################
#########################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(ComplexHeatmap)
  library(circlize)
  library(fs)
  library(cowplot)
  library(grid)
  library(pdftools)
  library(ggrepel)
})

# =============================================================================
# USER SETTINGS
# =============================================================================
SCATTER_PATHWAYS_TO_SHOW <- c(
  "Interferon alpha/beta signaling",
  "Interferon gamma signaling",
  "Interferon Signaling",
  "response to type II interferon",
  "response to interferon-beta",
  "cellular response to type II interferon",
  "Interferon gamma signaling",
  "response to type I interferon",
  "Interferon Signaling",
  "cellular response to type I interferon"
)

IN_FILE <- "PBMC_ILEUM_COMBINED_GSEA_ALL.csv"
OUT_DIR <- "PBMC_ILEUM_COMPARISON_PLOTS_FINE"
dir_create(OUT_DIR)

FDR_CUTOFF <- 0.05
TOP_N_HEATMAP <- 100
TOP_N_DELTA <- 20
TOP_N_CONCORDANT <- 100
TOP_N_SCATTER_LABELS <- 3
SCATTER_LABEL_WIDTH <- 40

GENESET_KEEP <- c("GO:BP", "REACTOME")

#CELLTYPES_KEEP <- c(
#  "B_cells", "CD4", "CD8", "Dendritic_cells",
#  "Monocytes", "NK_cells", "Progenitors", "T_cells"
#)

CELLTYPES_KEEP <- combined %>%
  dplyr::pull(cluster) %>%
  unique() %>%
  sort()

# =============================================================================
# HELPERS
# =============================================================================

clean_cluster_names <- function(x) {
  x %>%
    str_replace_all("\\+", "") %>%
    str_replace_all("__", "_") %>%
    str_replace_all(" ", "_") %>%
    str_replace_all("-", "_") %>%
    str_squish()
}

pretty_cluster_names <- function(x) {
  x %>%
    str_replace_all("_", " ")
}

save_plot <- function(plot_obj, filename, width, height, dpi = 300) {
  out <- file.path(OUT_DIR, filename)
  ggsave(out, plot_obj, width = width, height = height, dpi = dpi, limitsize = FALSE)
  out
}

make_pdf_panel <- function(file, label, title, dpi = 150) {
  img_raw <- pdftools::pdf_render_page(file, page = 1, dpi = dpi)
  
  grob_img <- tryCatch({
    rasterGrob(img_raw, interpolate = TRUE)
  }, error = function(e1) {
    tryCatch({
      rasterGrob(as.raster(img_raw), interpolate = TRUE)
    }, error = function(e2) {
      if (length(dim(img_raw)) == 3) {
        arr <- img_raw
        if (is.raw(arr)) arr <- as.integer(arr)
        if (max(arr, na.rm = TRUE) > 1) arr <- arr / 255
        rasterGrob(arr, interpolate = TRUE)
      } else {
        stop("Could not convert rendered PDF page into raster for: ", file)
      }
    })
  })
  
  ggdraw() +
    draw_label(
      label, x = 0.015, y = 0.985, hjust = 0, vjust = 1,
      fontface = "bold", size = 18
    ) +
    draw_label(
      title, x = 0.09, y = 0.98, hjust = 0, vjust = 1,
      fontface = "bold", size = 11
    ) +
    draw_grob(grob_img, x = 0.02, y = 0.02, width = 0.96, height = 0.90)
}

build_nes_matrix <- function(df_long, row_col = "row_id", col_col = "col_id", value_col = "NES", fill = 0) {
  mat_df <- df_long %>%
    dplyr::select(all_of(c(row_col, col_col, value_col))) %>%
    distinct() %>%
    pivot_wider(
      names_from = all_of(col_col),
      values_from = all_of(value_col),
      values_fill = fill
    ) %>%
    as.data.frame()
  
  rownames(mat_df) <- mat_df[[row_col]]
  mat_df[[row_col]] <- NULL
  as.matrix(mat_df)
}

make_dual_dataset_long <- function(df, include_direction = FALSE) {
  cols_keep <- c("pathway_id", "Description", "cluster", "geneset", "NES__PBMC", "NES__ILEUM")
  if (include_direction) cols_keep <- c(cols_keep, "direction_group")
  
  df %>%
    dplyr::select(all_of(cols_keep)) %>%
    pivot_longer(
      cols = c(NES__PBMC, NES__ILEUM),
      names_to = "dataset",
      values_to = "NES"
    ) %>%
    mutate(
      dataset = recode(dataset, NES__PBMC = "PBMC", NES__ILEUM = "DUODENUM"),
      col_id = paste(dataset, pretty_cluster_names(cluster), sep = " | "),
      row_id = paste0(Description, " [", geneset, "]")
    )
}

make_dotplot_long <- function(df, include_direction = FALSE) {
  cols_keep <- c("cluster", "Description", "NES__PBMC", "NES__ILEUM", "p.adjust__PBMC", "p.adjust__ILEUM")
  if (include_direction) cols_keep <- c(cols_keep, "direction_group")
  
  df %>%
    mutate(pathway_label = str_trunc(Description, 45)) %>%
    dplyr::select(all_of(cols_keep), pathway_label) %>%
    pivot_longer(
      cols = c(NES__PBMC, NES__ILEUM),
      names_to = "dataset",
      values_to = "NES"
    ) %>%
    mutate(
      dataset = recode(dataset, NES__PBMC = "PBMC", NES__ILEUM = "DUODENUM"),
      neglog10fdr = case_when(
        dataset == "PBMC" ~ -log10(p.adjust__PBMC),
        dataset == "DUODENUM" ~ -log10(p.adjust__ILEUM)
      )
    )
}

save_heatmap_pdf <- function(mat, filename, title, width, height,
                             row_split = NULL, column_split = NULL,
                             cluster_rows = TRUE, cluster_columns = TRUE,
                             row_font = 7, col_font = 8) {
  out <- file.path(OUT_DIR, filename)
  
  max_abs <- max(abs(mat), na.rm = TRUE)
  max_abs <- max(2, min(max_abs, 4))
  
  col_fun <- colorRamp2(
    c(-max_abs, 0, max_abs),
    c("#2166AC", "white", "#D6604D")
  )
  
  pdf(out, width = width, height = height)
  draw(
    Heatmap(
      mat,
      name             = "NES",
      col              = col_fun,
      cluster_rows     = cluster_rows,
      cluster_columns  = cluster_columns,
      row_split        = row_split,
      column_split     = column_split,
      row_names_gp     = gpar(fontsize = row_font),
      column_names_gp  = gpar(fontsize = col_font),
      column_title     = title,
      row_title        = "Pathways"
    )
  )
  dev.off()
  
  out
}

# =============================================================================
# LOAD
# =============================================================================

df <- fread(IN_FILE)

df <- df %>%
  mutate(
    cluster = clean_cluster_names(cluster),
    dataset = toupper(dataset),
    geneset = as.character(geneset),
    pathway_id = as.character(pathway_id),
    Description = as.character(Description),
    celltype_dataset = paste(dataset, cluster, sep = "_")
  ) %>% filter(geneset %in% GENESET_KEEP)

##FIX THIS
#if (!is.null(CELLTYPES_KEEP)) {
#  df <- df %>% filter(cluster %in% CELLTYPES_KEEP)
#}

df <- df %>% filter(dataset %in% c("PBMC", "ILEUM"))

# =============================================================================
# MATCH PBMC AND ILEUM COLUMN DATA, LABEL AS DUODENUM FOR OUTPUT
# =============================================================================

matched <- df %>%
  dplyr::select(
    pathway_id, ID, Description, cluster, geneset, dataset,
    NES, p.adjust
  ) %>%
  pivot_wider(
    names_from = dataset,
    values_from = c(NES, p.adjust),
    names_sep = "__"
  ) %>%
  filter(!is.na(NES__PBMC), !is.na(NES__ILEUM)) %>%
  mutate(
    delta_NES = NES__PBMC - NES__ILEUM,
    abs_delta_NES = abs(delta_NES),
    sig_status = case_when(
      p.adjust__PBMC <= FDR_CUTOFF & p.adjust__ILEUM <= FDR_CUTOFF ~ "Both significant",
      p.adjust__PBMC <= FDR_CUTOFF & p.adjust__ILEUM > FDR_CUTOFF  ~ "PBMC only",
      p.adjust__PBMC > FDR_CUTOFF  & p.adjust__ILEUM <= FDR_CUTOFF ~ "DUODENUM only",
      TRUE ~ "Neither significant"
    ),
    direction_group = case_when(
      NES__PBMC > 0 & NES__ILEUM > 0 ~ "Both Up",
      NES__PBMC < 0 & NES__ILEUM < 0 ~ "Both Down",
      NES__PBMC > 0 & NES__ILEUM < 0 ~ "Opposite: PBMC Up / Duodenum Down",
      NES__PBMC < 0 & NES__ILEUM > 0 ~ "Opposite: PBMC Down / Duodenum Up",
      TRUE ~ "Mixed/Zero"
    ),
    quadrant = case_when(
      NES__PBMC > 0 & NES__ILEUM > 0 ~ "Q1: Both Up",
      NES__PBMC < 0 & NES__ILEUM > 0 ~ "Q2: Duodenum Up / PBMC Down",
      NES__PBMC < 0 & NES__ILEUM < 0 ~ "Q3: Both Down",
      NES__PBMC > 0 & NES__ILEUM < 0 ~ "Q4: PBMC Up / Duodenum Down",
      TRUE ~ "On Axis"
    ),
    concordance_strength = abs(NES__PBMC) + abs(NES__ILEUM),
    mean_NES = (NES__PBMC + NES__ILEUM) / 2
  )

fwrite(matched, file.path(OUT_DIR, "matched_pbmc_duodenum_pathways.csv"))

fwrite(
  matched %>%
    dplyr::select(
      pathway_id, ID, Description, cluster, geneset,
      NES__PBMC, NES__ILEUM, p.adjust__PBMC, p.adjust__ILEUM,
      sig_status, quadrant, delta_NES, abs_delta_NES
    ),
  file.path(OUT_DIR, "scatter_points_with_quadrants.csv")
)

quad_counts <- matched %>%
  count(cluster, quadrant, sort = TRUE)

fwrite(quad_counts, file.path(OUT_DIR, "quadrant_counts.csv"))

# =============================================================================
# TOP CONCORDANT
# =============================================================================

top_both_up <- matched %>%
  filter(direction_group == "Both Up") %>%
  arrange(desc(concordance_strength), desc(mean_NES)) %>%
  slice_head(n = TOP_N_CONCORDANT)

top_both_down <- matched %>%
  filter(direction_group == "Both Down") %>%
  arrange(desc(concordance_strength), mean_NES) %>%
  slice_head(n = TOP_N_CONCORDANT)

top_concordant <- bind_rows(top_both_up, top_both_down)

fwrite(top_both_up, file.path(OUT_DIR, "top100_both_up_pathways.csv"))
fwrite(top_both_down, file.path(OUT_DIR, "top100_both_down_pathways.csv"))
fwrite(top_concordant, file.path(OUT_DIR, "top100_both_up_and_down_combined.csv"))

# =============================================================================
# HEATMAP 1
# =============================================================================

heatmap_df <- matched %>%
  group_by(pathway_id, Description, geneset) %>%
  summarise(mean_abs_delta = mean(abs_delta_NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_delta)) %>%
  slice_head(n = TOP_N_HEATMAP)

heatmap_long <- matched %>%
  filter(pathway_id %in% heatmap_df$pathway_id) %>%
  make_dual_dataset_long()

heatmap_mat <- build_nes_matrix(heatmap_long)

heatmap_top_variable_pdf <- save_heatmap_pdf(
  mat = heatmap_mat,
  filename = "heatmap_top_variable_shared_pathways.pdf",
  title = "PBMC vs Duodenum shared pathways",
  width = 14,
  height = 16
)

# =============================================================================
# HEATMAP 2
# =============================================================================

concordant_heatmap_long <- make_dual_dataset_long(top_concordant, include_direction = TRUE)
concordant_heatmap_mat <- build_nes_matrix(concordant_heatmap_long)

row_split_df <- top_concordant %>%
  distinct(pathway_id, Description, geneset, direction_group) %>%
  mutate(row_id = paste0(Description, " [", geneset, "]")) %>%
  distinct(row_id, direction_group)

row_split_vec <- row_split_df$direction_group[
  match(rownames(concordant_heatmap_mat), row_split_df$row_id)
]
row_split_vec <- factor(row_split_vec, levels = c("Both Up", "Both Down"))

heatmap_concordant_pdf <- save_heatmap_pdf(
  mat = concordant_heatmap_mat,
  filename = "heatmap_top100_both_up_and_down.pdf",
  title = "Top 100 pathways both up and top 100 pathways both down",
  width = 14,
  height = 18,
  row_split = row_split_vec
)

# =============================================================================
# HEATMAP 3 — Top 100 discordant pathways (opposite directions)
# =============================================================================
top_discordant <- matched %>%
  filter(
    sig_status != "Neither significant",
    quadrant %in% c(
      "Q2: Duodenum Up / PBMC Down",
      "Q4: PBMC Up / Duodenum Down"
    )
  ) %>%
  mutate(
    direction_group = case_when(
      quadrant == "Q4: PBMC Up / Duodenum Down" ~ "PBMC Up / Duodenum Down",
      quadrant == "Q2: Duodenum Up / PBMC Down" ~ "Duodenum Up / PBMC Down"
    ),
    # rank by how extreme the discordance is — largest absolute NES in each direction
    discordance_score = abs(NES__PBMC) + abs(NES__ILEUM)
  ) %>%
  group_by(direction_group, geneset) %>%
  slice_max(order_by = discordance_score, n = 50, with_ties = FALSE) %>%
  ungroup()

discordant_heatmap_long <- make_dual_dataset_long(top_discordant, include_direction = TRUE)
discordant_heatmap_mat  <- build_nes_matrix(discordant_heatmap_long)

row_split_df_disc <- top_discordant %>%
  distinct(pathway_id, Description, geneset, direction_group) %>%
  mutate(row_id = paste0(Description, " [", geneset, "]")) %>%
  distinct(row_id, direction_group)

row_split_vec_disc <- row_split_df_disc$direction_group[
  match(rownames(discordant_heatmap_mat), row_split_df_disc$row_id)
]
row_split_vec_disc <- factor(
  row_split_vec_disc,
  levels = c("PBMC Up / Duodenum Down", "Duodenum Up / PBMC Down")
)

heatmap_discordant_paired_pdf <- save_heatmap_pdf(
  mat              = discordant_heatmap_mat_paired,
  filename         = "heatmap_top100_discordant_paired.pdf",
  title            = "Top 100 pathways PBMC Up/Duodenum Down + Top 100 Duodenum Up/PBMC Down\n(cell types paired by column)",
  width            = 14,
  height           = 18,
  row_split        = row_split_vec_disc,
  column_split     = col_split_vec,
  cluster_columns  = FALSE
)


# =============================================================================
# DELTA NES BARPLOT
# =============================================================================

top_delta_by_cluster <- matched %>%
  group_by(cluster) %>%
  slice_max(order_by = abs_delta_NES, n = TOP_N_DELTA, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    pathway_label = str_trunc(Description, 50),
    pathway_label = fct_reorder(pathway_label, delta_NES)
  )

p_delta <- ggplot(
  top_delta_by_cluster,
  aes(x = delta_NES, y = pathway_label, fill = sig_status)
) +
  geom_col() +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    title = "Top PBMC vs Duodenum pathway differences",
    subtitle = "delta_NES = PBMC NES - Duodenum NES",
    x = "delta NES",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

delta_pdf <- save_plot(p_delta, "deltaNES_top_pathways_by_cluster.pdf", 16, 12)

# =============================================================================
# SCATTER PLOT WITH BIGGER SIZE + LESS OVERLAP + USER LIST + AUTO TOP
# =============================================================================

scatter_df <- matched %>%
  filter(sig_status != "Neither significant")

top_points_auto <- scatter_df %>%
  mutate(label_strength = abs(NES__PBMC) + abs(NES__ILEUM)) %>%
  group_by(cluster, quadrant) %>%
  slice_max(order_by = label_strength, n = TOP_N_SCATTER_LABELS, with_ties = FALSE) %>%
  ungroup()

if (!is.null(SCATTER_PATHWAYS_TO_SHOW)) {
  top_points_user <- scatter_df %>%
    filter(
      pathway_id %in% SCATTER_PATHWAYS_TO_SHOW |
        Description %in% SCATTER_PATHWAYS_TO_SHOW
    )
} else {
  top_points_user <- scatter_df %>% slice(0)
}

top_points <- bind_rows(top_points_auto, top_points_user) %>%
  distinct(cluster, pathway_id, .keep_all = TRUE)

quadrant_colors <- c(
  "Q1: Both Up" = "#D6604D",
  "Q2: Duodenum Up / PBMC Down" = "#4393C3",
  "Q3: Both Down" = "#2166AC",
  "Q4: PBMC Up / Duodenum Down" = "#B2182B",
  "On Axis" = "grey70"
)

sig_shapes <- c(
  "Both significant" = 16,
  "PBMC only" = 17,
  "DUODENUM only" = 15
)

p_scatter <- scatter_df %>%
  ggplot(aes(
    x = NES__PBMC,
    y = NES__ILEUM,
    color = quadrant,
    shape = sig_status
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(
    data = top_points,
    aes(label = str_trunc(Description, SCATTER_LABEL_WIDTH)),
    size = 2.8,
    max.overlaps = 100,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.2,
    force = 2,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  facet_wrap(~ cluster, scales = "free", ncol = 3) +
  scale_color_manual(values = quadrant_colors) +
  scale_shape_manual(values = sig_shapes) +
  labs(
    title = "PBMC vs Duodenum NES correlation with quadrant classification",
    x = "PBMC NES",
    y = "Duodenum NES",
    color = "Quadrant",
    shape = "Significance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

scatter_pdf <- save_plot(
  p_scatter,
  "scatter_PBMC_vs_Duodenum_NES.pdf",
  26,
  20
)

save_plot(
  p_scatter,
  "scatter_PBMC_vs_Duodenum_NES.png",
  26,
  20
)

p_scatter_disagreement <- matched %>%
  filter(sig_status != "Neither significant") %>%
  filter(quadrant %in% c("Q2: Duodenum Up / PBMC Down", "Q4: PBMC Up / Duodenum Down")) %>%
  ggplot(aes(x = NES__PBMC, y = NES__ILEUM, color = quadrant, shape = sig_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  geom_point(alpha = 0.7, size = 1.4) +
  facet_wrap(~ cluster, scales = "free", ncol = 3) +
  scale_color_manual(values = quadrant_colors) +
  scale_shape_manual(values = sig_shapes) +
  labs(
    title = "PBMC vs Duodenum disagreement pathways",
    x = "PBMC NES",
    y = "Duodenum NES",
    color = "Quadrant",
    shape = "Significance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

save_plot(p_scatter_disagreement, "scatter_PBMC_vs_Duodenum_disagreement_only.pdf", 18, 12)

# =============================================================================
# DOT PLOT: SIGNIFICANT SHARED PATHWAYS
# =============================================================================

dot_df <- matched %>%
  filter(sig_status != "Neither significant") %>%
  group_by(cluster) %>%
  slice_max(order_by = abs_delta_NES, n = TOP_N_DELTA, with_ties = FALSE) %>%
  ungroup() %>%
  make_dotplot_long()

p_dot <- ggplot(
  dot_df,
  aes(x = dataset, y = pathway_label, size = neglog10fdr, color = NES)
) +
  geom_point() +
  facet_wrap(~ cluster, scales = "free_y") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0) +
  labs(
    title = "PBMC vs Duodenum significant pathway comparison",
    x = NULL,
    y = NULL,
    size = "-log10(FDR)",
    color = "NES"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dot_sig_pdf <- save_plot(p_dot, "dotplot_significant_pathways_PBMC_vs_Duodenum.pdf", 16, 12)

# =============================================================================
# DOT PLOT: TOP CONCORDANT
# =============================================================================

dot_df_concordant <- make_dotplot_long(top_concordant, include_direction = TRUE)

p_dot_concordant <- ggplot(
  dot_df_concordant,
  aes(x = dataset, y = pathway_label, size = neglog10fdr, color = NES)
) +
  geom_point() +
  facet_grid(direction_group ~ cluster, scales = "free_y", space = "free_y") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0) +
  labs(
    title = "Top 100 concordant up and top 100 concordant down pathways",
    x = NULL,
    y = NULL,
    size = "-log10(FDR)",
    color = "NES"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dot_concordant_pdf <- save_plot(p_dot_concordant, "dotplot_top100_both_up_and_down.pdf", 18, 14)

# =============================================================================
# SUMMARY PLOTS
# =============================================================================

summary_counts <- matched %>%
  group_by(cluster, sig_status) %>%
  summarise(n = n(), .groups = "drop")

p_counts <- ggplot(summary_counts, aes(x = cluster, y = n, fill = sig_status)) +
  geom_col(position = "stack") +
  labs(
    title = "Number of matched pathways by significance pattern",
    x = "Cell type",
    y = "Number of pathways"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_sig_pdf <- save_plot(p_counts, "summary_counts_sigstatus.pdf", 10, 6)

summary_concordant <- matched %>%
  filter(direction_group %in% c("Both Up", "Both Down")) %>%
  group_by(cluster, direction_group) %>%
  summarise(n = n(), .groups = "drop")

p_concordant_counts <- ggplot(summary_concordant, aes(x = cluster, y = n, fill = direction_group)) +
  geom_col(position = "stack") +
  labs(
    title = "Number of concordant pathways by cell type",
    x = "Cell type",
    y = "Number of pathways"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_concordant_pdf <- save_plot(p_concordant_counts, "summary_counts_both_up_both_down.pdf", 10, 6)

p_quadrant_counts <- ggplot(quad_counts, aes(x = cluster, y = n, fill = quadrant)) +
  geom_col(position = "stack") +
  labs(
    title = "Number of pathways by scatter quadrant",
    x = "Cell type",
    y = "Number of pathways"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_quadrant_counts, "summary_counts_quadrants.pdf", 11, 6)

# =============================================================================
# FINAL PANEL
# =============================================================================

panel_info <- list(
  list(file = scatter_pdf, label = "A", title = "PBMC vs Duodenum NES correlation"),
  list(file = delta_pdf, label = "B", title = "Top PBMC vs Duodenum pathway differences"),
  list(file = dot_concordant_pdf, label = "C", title = "Top 100 concordant up/down pathways"),
  list(file = heatmap_concordant_pdf, label = "D", title = "Concordant pathway heatmap"),
  list(file = summary_sig_pdf, label = "E", title = "Matched pathways by significance pattern"),
  list(file = summary_concordant_pdf, label = "F", title = "Concordant pathways by cell type")
)

panels <- lapply(panel_info, function(x) {
  make_pdf_panel(x$file, x$label, x$title)
})

title_panel <- ggdraw() +
  draw_label(
    "PBMC versus Duodenum immune-cell pathway comparison",
    x = 0.5, y = 0.70,
    hjust = 0.5, vjust = 0.5,
    fontface = "bold", size = 18
  ) +
  draw_label(
    "Integrated multi-panel summary of concordance, difference, and significance across matched Monaco cell types",
    x = 0.5, y = 0.25,
    hjust = 0.5, vjust = 0.5,
    size = 11
  )

top_row <- plot_grid(panels[[1]], panels[[2]], ncol = 2, align = "hv")
mid_row <- plot_grid(panels[[3]], panels[[4]], ncol = 2, align = "hv")
bot_row <- plot_grid(panels[[5]], panels[[6]], ncol = 2, align = "hv")

final_panel <- plot_grid(
  title_panel,
  top_row,
  mid_row,
  bot_row,
  ncol = 1,
  rel_heights = c(0.10, 1, 1, 1)
)

save_plot(final_panel, "PBMC_DUODENUM_paper_style_panel.pdf", 16, 22)
save_plot(final_panel, "PBMC_DUODENUM_paper_style_panel.png", 16, 22)

cat("Done. Plots written to:", OUT_DIR, "\n")

#####################



suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(ComplexHeatmap)
  library(circlize)
  library(fs)
  library(cowplot)
  library(grid)
  library(pdftools)
  library(ggrepel)
})

# =============================================================================
# USER SETTINGS
# =============================================================================
SCATTER_PATHWAYS_TO_SHOW <- c("Interferon alpha/beta signaling", "Interferon gamma signaling", "Interferon Signaling", "response to type II interferon", "response to interferon-beta", "cellular response to type II interferon", "Interferon gamma signaling", "response to type I interferon", "Interferon Signaling", "cellular response to type I interferon")

IN_FILE <- "PBMC_ILEUM_COMBINED_GSEA_ALL_FINE_ANN.csv"
OUT_DIR <- "PBMC_ILEUM_COMPARISON_PLOTS_FINE"
dir_create(OUT_DIR)


FDR_CUTOFF <- 0.05
TOP_N_HEATMAP <- 100
TOP_N_DELTA <- 20
TOP_N_CONCORDANT <- 100
TOP_N_SCATTER_LABELS <- 3
SCATTER_LABEL_WIDTH <- 40

GENESET_KEEP <- c("GO:BP", "REACTOME")


CELLTYPES_KEEP <- combined %>%
  dplyr::pull(cluster) %>%
  unique() %>%
  sort()

# =============================================================================
# HELPERS
# =============================================================================

clean_cluster_names <- function(x) {
  x %>%
    str_replace_all("\\+", "") %>%
    str_replace_all("__", "_") %>%
    str_replace_all(" ", "_") %>%
    str_replace_all("-", "_") %>%
    str_squish()
}

pretty_cluster_names <- function(x) {
  x %>%
    str_replace_all("_", " ")
}

save_plot <- function(plot_obj, filename, width, height, dpi = 300) {
  out <- file.path(OUT_DIR, filename)
  ggsave(out, plot_obj, width = width, height = height, dpi = dpi, limitsize = FALSE)
  out
}

make_pdf_panel <- function(file, label, title, dpi = 150) {
  img_raw <- pdftools::pdf_render_page(file, page = 1, dpi = dpi)
  
  grob_img <- tryCatch({
    rasterGrob(img_raw, interpolate = TRUE)
  }, error = function(e1) {
    tryCatch({
      rasterGrob(as.raster(img_raw), interpolate = TRUE)
    }, error = function(e2) {
      if (length(dim(img_raw)) == 3) {
        arr <- img_raw
        if (is.raw(arr)) arr <- as.integer(arr)
        if (max(arr, na.rm = TRUE) > 1) arr <- arr / 255
        rasterGrob(arr, interpolate = TRUE)
      } else {
        stop("Could not convert rendered PDF page into raster for: ", file)
      }
    })
  })
  
  ggdraw() +
    draw_label(label, x = 0.015, y = 0.985, hjust = 0, vjust = 1,
               fontface = "bold", size = 18) +
    draw_label(title, x = 0.09, y = 0.98, hjust = 0, vjust = 1,
               fontface = "bold", size = 11) +
    draw_grob(grob_img, x = 0.02, y = 0.02, width = 0.96, height = 0.90)
}

build_nes_matrix <- function(df_long, row_col = "row_id", col_col = "col_id", value_col = "NES", fill = 0) {
  mat_df <- df_long %>%
    dplyr::select(all_of(c(row_col, col_col, value_col))) %>%
    distinct() %>%
    pivot_wider(
      names_from = all_of(col_col),
      values_from = all_of(value_col),
      values_fill = fill
    ) %>%
    as.data.frame()
  
  rownames(mat_df) <- mat_df[[row_col]]
  mat_df[[row_col]] <- NULL
  as.matrix(mat_df)
}

make_dual_dataset_long <- function(df, include_direction = FALSE) {
  cols_keep <- c("pathway_id", "Description", "cluster", "geneset", "NES__PBMC", "NES__ILEUM")
  if (include_direction) cols_keep <- c(cols_keep, "direction_group")
  
  df %>%
    dplyr::select(all_of(cols_keep)) %>%
    pivot_longer(
      cols = c(NES__PBMC, NES__ILEUM),
      names_to = "dataset",
      values_to = "NES"
    ) %>%
    mutate(
      dataset = recode(dataset, NES__PBMC = "PBMC", NES__ILEUM = "ILEUM"),
      col_id = paste(dataset, pretty_cluster_names(cluster), sep = " | "),
      row_id = paste0(Description, " [", geneset, "]")
    )
}

make_dotplot_long <- function(df, include_direction = FALSE) {
  cols_keep <- c("cluster", "Description", "NES__PBMC", "NES__ILEUM", "p.adjust__PBMC", "p.adjust__ILEUM")
  if (include_direction) cols_keep <- c(cols_keep, "direction_group")
  
  df %>%
    mutate(pathway_label = str_trunc(Description, 45)) %>%
    dplyr::select(all_of(cols_keep), pathway_label) %>%
    pivot_longer(
      cols = c(NES__PBMC, NES__ILEUM),
      names_to = "dataset",
      values_to = "NES"
    ) %>%
    mutate(
      dataset = recode(dataset, NES__PBMC = "PBMC", NES__ILEUM = "ILEUM"),
      neglog10fdr = case_when(
        dataset == "PBMC" ~ -log10(p.adjust__PBMC),
        dataset == "ILEUM" ~ -log10(p.adjust__ILEUM)
      )
    )
}

save_heatmap_pdf <- function(mat, filename, title, width, height,
                             row_split = NULL, cluster_rows = TRUE, cluster_columns = TRUE,
                             row_font = 7, col_font = 8) {
  out <- file.path(OUT_DIR, filename)
  
  max_abs <- max(abs(mat), na.rm = TRUE)
  max_abs <- max(2, min(max_abs, 4))
  
  col_fun <- colorRamp2(
    c(-max_abs, 0, max_abs),
    c("#2166AC", "white", "#D6604D")
  )
  
  pdf(out, width = width, height = height)
  draw(
    Heatmap(
      mat,
      name = "NES",
      col = col_fun,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      row_split = row_split,
      row_names_gp = gpar(fontsize = row_font),
      column_names_gp = gpar(fontsize = col_font),
      column_title = title,
      row_title = "Pathways"
    )
  )
  dev.off()
  
  out
}

# =============================================================================
# LOAD
# =============================================================================

df <- fread(IN_FILE)

df <- df %>%
  mutate(
    cluster = clean_cluster_names(cluster),
    dataset = toupper(dataset),
    geneset = as.character(geneset),
    pathway_id = as.character(pathway_id),
    Description = as.character(Description),
    celltype_dataset = paste(dataset, cluster, sep = "_")
  ) %>%
  filter(geneset %in% GENESET_KEEP)

#if (!is.null(CELLTYPES_KEEP)) {
#  df <- df %>% filter(cluster %in% CELLTYPES_KEEP)
#}

df <- df %>% filter(dataset %in% c("PBMC", "ILEUM"))

# =============================================================================
# MATCH PBMC AND ILEUM
# =============================================================================

matched <- df %>%
  dplyr::select(
    pathway_id, ID, Description, cluster, geneset, dataset,
    NES, p.adjust
  ) %>%
  pivot_wider(
    names_from = dataset,
    values_from = c(NES, p.adjust),
    names_sep = "__"
  ) %>%
  filter(!is.na(NES__PBMC), !is.na(NES__ILEUM)) %>%
  mutate(
    delta_NES = NES__PBMC - NES__ILEUM,
    abs_delta_NES = abs(delta_NES),
    sig_status = case_when(
      p.adjust__PBMC <= FDR_CUTOFF & p.adjust__ILEUM <= FDR_CUTOFF ~ "Both significant",
      p.adjust__PBMC <= FDR_CUTOFF & p.adjust__ILEUM > FDR_CUTOFF  ~ "PBMC only",
      p.adjust__PBMC > FDR_CUTOFF  & p.adjust__ILEUM <= FDR_CUTOFF ~ "ILEUM only",
      TRUE ~ "Neither significant"
    ),
    direction_group = case_when(
      NES__PBMC > 0 & NES__ILEUM > 0 ~ "Both Up",
      NES__PBMC < 0 & NES__ILEUM < 0 ~ "Both Down",
      NES__PBMC > 0 & NES__ILEUM < 0 ~ "Opposite: PBMC Up / Ileum Down",
      NES__PBMC < 0 & NES__ILEUM > 0 ~ "Opposite: PBMC Down / Ileum Up",
      TRUE ~ "Mixed/Zero"
    ),
    quadrant = case_when(
      NES__PBMC > 0 & NES__ILEUM > 0 ~ "Q1: Both Up",
      NES__PBMC < 0 & NES__ILEUM > 0 ~ "Q2: Ileum Up / PBMC Down",
      NES__PBMC < 0 & NES__ILEUM < 0 ~ "Q3: Both Down",
      NES__PBMC > 0 & NES__ILEUM < 0 ~ "Q4: PBMC Up / Ileum Down",
      TRUE ~ "On Axis"
    ),
    concordance_strength = abs(NES__PBMC) + abs(NES__ILEUM),
    mean_NES = (NES__PBMC + NES__ILEUM) / 2
  )

fwrite(matched, file.path(OUT_DIR, "matched_pbmc_ileum_pathways.csv"))

fwrite(
  matched %>%
    dplyr::select(pathway_id, ID, Description, cluster, geneset,
           NES__PBMC, NES__ILEUM, p.adjust__PBMC, p.adjust__ILEUM,
           sig_status, quadrant, delta_NES, abs_delta_NES),
  file.path(OUT_DIR, "scatter_points_with_quadrants.csv")
)

quad_counts <- matched %>%
  count(cluster, quadrant, sort = TRUE)

fwrite(quad_counts, file.path(OUT_DIR, "quadrant_counts.csv"))

# =============================================================================
# TOP CONCORDANT
# =============================================================================

top_both_up <- matched %>%
  filter(direction_group == "Both Up") %>%
  arrange(desc(concordance_strength), desc(mean_NES)) %>%
  slice_head(n = TOP_N_CONCORDANT)

top_both_down <- matched %>%
  filter(direction_group == "Both Down") %>%
  arrange(desc(concordance_strength), mean_NES) %>%
  slice_head(n = TOP_N_CONCORDANT)

top_concordant <- bind_rows(top_both_up, top_both_down)

fwrite(top_both_up, file.path(OUT_DIR, "top100_both_up_pathways.csv"))
fwrite(top_both_down, file.path(OUT_DIR, "top100_both_down_pathways.csv"))
fwrite(top_concordant, file.path(OUT_DIR, "top100_both_up_and_down_combined.csv"))

# =============================================================================
# HEATMAP 1
# =============================================================================

heatmap_df <- matched %>%
  group_by(pathway_id, Description, geneset) %>%
  summarise(mean_abs_delta = mean(abs_delta_NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_delta)) %>%
  slice_head(n = TOP_N_HEATMAP)

heatmap_long <- matched %>%
  dplyr::filter(pathway_id %in% heatmap_df$pathway_id) %>%
  make_dual_dataset_long()

heatmap_mat <- build_nes_matrix(heatmap_long)

heatmap_top_variable_pdf <- save_heatmap_pdf(
  mat = heatmap_mat,
  filename = "heatmap_top_variable_shared_pathways.pdf",
  title = "PBMC vs Ileum shared pathways",
  width = 14,
  height = 16
)

# =============================================================================
# HEATMAP 2
# =============================================================================

concordant_heatmap_long <- make_dual_dataset_long(top_concordant, include_direction = TRUE)
concordant_heatmap_mat <- build_nes_matrix(concordant_heatmap_long)

row_split_df <- top_concordant %>%
  distinct(pathway_id, Description, geneset, direction_group) %>%
  mutate(row_id = paste0(Description, " [", geneset, "]")) %>%
  distinct(row_id, direction_group)

row_split_vec <- row_split_df$direction_group[
  match(rownames(concordant_heatmap_mat), row_split_df$row_id)
]
row_split_vec <- factor(row_split_vec, levels = c("Both Up", "Both Down"))

heatmap_concordant_pdf <- save_heatmap_pdf(
  mat = concordant_heatmap_mat,
  filename = "heatmap_top100_both_up_and_down.pdf",
  title = "Top 100 pathways both up and top 100 pathways both down",
  width = 14,
  height = 18,
  row_split = row_split_vec
)


# =============================================================================
# HEATMAP 3 — Top 100 discordant pathways (opposite directions)
# =============================================================================
top_discordant <- matched %>%
  filter(
    sig_status != "Neither significant",
    quadrant %in% c(
      "Q2: Duodenum Up / PBMC Down",
      "Q4: PBMC Up / Duodenum Down"
    )
  ) %>%
  mutate(
    direction_group = case_when(
      quadrant == "Q4: PBMC Up / Duodenum Down" ~ "PBMC Up / Duodenum Down",
      quadrant == "Q2: Duodenum Up / PBMC Down" ~ "Duodenum Up / PBMC Down"
    ),
    # rank by how extreme the discordance is — largest absolute NES in each direction
    discordance_score = abs(NES__PBMC) + abs(NES__ILEUM)
  ) %>%
  group_by(direction_group, geneset) %>%
  slice_max(order_by = discordance_score, n = 50, with_ties = FALSE) %>%
  ungroup()

discordant_heatmap_long <- make_dual_dataset_long(top_discordant, include_direction = TRUE)
discordant_heatmap_mat  <- build_nes_matrix(discordant_heatmap_long)

row_split_df_disc <- top_discordant %>%
  distinct(pathway_id, Description, geneset, direction_group) %>%
  mutate(row_id = paste0(Description, " [", geneset, "]")) %>%
  distinct(row_id, direction_group)

row_split_vec_disc <- row_split_df_disc$direction_group[
  match(rownames(discordant_heatmap_mat), row_split_df_disc$row_id)
]
row_split_vec_disc <- factor(
  row_split_vec_disc,
  levels = c("PBMC Up / Duodenum Down", "Duodenum Up / PBMC Down")
)

heatmap_discordant_pdf <- save_heatmap_pdf(
  mat        = discordant_heatmap_mat,
  filename   = "heatmap_top100_discordant.pdf",
  title      = "Top 50 pathways PBMC Up/Duodenum Down + Top 50 pathways Duodenum Up/PBMC Down",
  width      = 14,
  height     = 18,
  row_split  = row_split_vec_disc
)

# =============================================================================
# DELTA NES BARPLOT
# =============================================================================

top_delta_by_cluster <- matched %>%
  group_by(cluster) %>%
  slice_max(order_by = abs_delta_NES, n = TOP_N_DELTA, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    pathway_label = str_trunc(Description, 50),
    pathway_label = fct_reorder(pathway_label, delta_NES)
  )

p_delta <- ggplot(top_delta_by_cluster,
                  aes(x = delta_NES, y = pathway_label, fill = sig_status)) +
  geom_col() +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    title = "Top PBMC vs Ileum pathway differences",
    subtitle = "delta_NES = PBMC NES - Ileum NES",
    x = "delta NES",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

delta_pdf <- save_plot(p_delta, "deltaNES_top_pathways_by_cluster.pdf", 16, 12)

# =============================================================================
# SCATTER PLOT WITH BIGGER SIZE + LESS OVERLAP
# =============================================================================

# =============================================================================
# SCATTER PLOT WITH BIGGER SIZE + LESS OVERLAP + USER LIST + AUTO TOP
# =============================================================================

scatter_df <- matched %>%
  filter(sig_status != "Neither significant")

# automatic top labels
top_points_auto <- scatter_df %>%
  mutate(label_strength = abs(NES__PBMC) + abs(NES__ILEUM)) %>%
  group_by(cluster, quadrant) %>%
  slice_max(order_by = label_strength, n = TOP_N_SCATTER_LABELS, with_ties = FALSE) %>%
  ungroup()

# user-specified labels
if (!is.null(SCATTER_PATHWAYS_TO_SHOW)) {
  top_points_user <- scatter_df %>%
    filter(
      pathway_id %in% SCATTER_PATHWAYS_TO_SHOW |
        Description %in% SCATTER_PATHWAYS_TO_SHOW
    )
} else {
  top_points_user <- scatter_df %>% slice(0)
}

# combine both and remove duplicates
top_points <- bind_rows(top_points_auto, top_points_user) %>%
  distinct(cluster, pathway_id, .keep_all = TRUE)

quadrant_colors <- c(
  "Q1: Both Up" = "#D6604D",
  "Q2: Ileum Up / PBMC Down" = "#4393C3",
  "Q3: Both Down" = "#2166AC",
  "Q4: PBMC Up / Ileum Down" = "#B2182B",
  "On Axis" = "grey70"
)

sig_shapes <- c(
  "Both significant" = 16,
  "PBMC only" = 17,
  "ILEUM only" = 15
)

p_scatter <- scatter_df %>%
  ggplot(aes(
    x = NES__PBMC,
    y = NES__ILEUM,
    color = quadrant,
    shape = sig_status
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(
    data = top_points,
    aes(label = str_trunc(Description, SCATTER_LABEL_WIDTH)),
    size = 2.8,
    max.overlaps = 100,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.2,
    force = 2,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  facet_wrap(~ cluster, scales = "free", ncol = 3) +
  scale_color_manual(values = quadrant_colors) +
  scale_shape_manual(values = sig_shapes) +
  labs(
    title = "PBMC vs Duodenum NES correlation with quadrant classification",
    x = "PBMC NES",
    y = "Duodenum NES",
    color = "Quadrant",
    shape = "Significance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

scatter_pdf <- save_plot(
  p_scatter,
  "scatter_PBMC_vs_Duodenum_NES.pdf",
  26,
  20
)

save_plot(
  p_scatter,
  "scatter_PBMC_vs_Duodenum_NES.png",
  26,
  20
)

p_scatter_disagreement <- matched %>%
  filter(sig_status != "Neither significant") %>%
  filter(quadrant %in% c("Q2: Ileum Up / PBMC Down", "Q4: PBMC Up / Ileum Down")) %>%
  ggplot(aes(x = NES__PBMC, y = NES__ILEUM, color = quadrant, shape = sig_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  geom_point(alpha = 0.7, size = 1.4) +
  facet_wrap(~ cluster, scales = "free", ncol = 3) +
  scale_color_manual(values = quadrant_colors) +
  scale_shape_manual(values = sig_shapes) +
  labs(
    title = "PBMC vs Ileum disagreement pathways",
    x = "PBMC NES",
    y = "Duodenum NES",
    color = "Quadrant",
    shape = "Significance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

save_plot(p_scatter_disagreement, "scatter_PBMC_vs_ILEUM_disagreement_only.pdf", 18, 12)

# =============================================================================
# DOT PLOT: SIGNIFICANT SHARED PATHWAYS
# =============================================================================

dot_df <- matched %>%
  filter(sig_status != "Neither significant") %>%
  group_by(cluster) %>%
  slice_max(order_by = abs_delta_NES, n = TOP_N_DELTA, with_ties = FALSE) %>%
  ungroup() %>%
  make_dotplot_long()

p_dot <- ggplot(dot_df,
                aes(x = dataset, y = pathway_label, size = neglog10fdr, color = NES)) +
  geom_point() +
  facet_wrap(~ cluster, scales = "free_y") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0) +
  labs(
    title = "PBMC vs Ileum significant pathway comparison",
    x = NULL,
    y = NULL,
    size = "-log10(FDR)",
    color = "NES"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dot_sig_pdf <- save_plot(p_dot, "dotplot_significant_pathways_PBMC_vs_ILEUM.pdf", 16, 12)

# =============================================================================
# DOT PLOT: TOP CONCORDANT
# =============================================================================

dot_df_concordant <- make_dotplot_long(top_concordant, include_direction = TRUE)

p_dot_concordant <- ggplot(
  dot_df_concordant,
  aes(x = dataset, y = pathway_label, size = neglog10fdr, color = NES)
) +
  geom_point() +
  facet_grid(direction_group ~ cluster, scales = "free_y", space = "free_y") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0) +
  labs(
    title = "Top 100 concordant up and top 100 concordant down pathways",
    x = NULL,
    y = NULL,
    size = "-log10(FDR)",
    color = "NES"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dot_concordant_pdf <- save_plot(p_dot_concordant, "dotplot_top100_both_up_and_down.pdf", 18, 14)

# =============================================================================
# SUMMARY PLOTS
# =============================================================================

summary_counts <- matched %>%
  group_by(cluster, sig_status) %>%
  summarise(n = n(), .groups = "drop")

p_counts <- ggplot(summary_counts,
                   aes(x = cluster, y = n, fill = sig_status)) +
  geom_col(position = "stack") +
  labs(
    title = "Number of matched pathways by significance pattern",
    x = "Cell type",
    y = "Number of pathways"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_sig_pdf <- save_plot(p_counts, "summary_counts_sigstatus.pdf", 10, 6)

summary_concordant <- matched %>%
  filter(direction_group %in% c("Both Up", "Both Down")) %>%
  group_by(cluster, direction_group) %>%
  summarise(n = n(), .groups = "drop")

p_concordant_counts <- ggplot(summary_concordant,
                              aes(x = cluster, y = n, fill = direction_group)) +
  geom_col(position = "stack") +
  labs(
    title = "Number of concordant pathways by cell type",
    x = "Cell type",
    y = "Number of pathways"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_concordant_pdf <- save_plot(p_concordant_counts, "summary_counts_both_up_both_down.pdf", 10, 6)

p_quadrant_counts <- ggplot(quad_counts, aes(x = cluster, y = n, fill = quadrant)) +
  geom_col(position = "stack") +
  labs(
    title = "Number of pathways by scatter quadrant",
    x = "Cell type",
    y = "Number of pathways"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_quadrant_counts, "summary_counts_quadrants.pdf", 11, 6)

# =============================================================================
# FINAL PANEL
# =============================================================================

panel_info <- list(
  list(file = scatter_pdf, label = "A", title = "PBMC vs Ileum NES correlation"),
  list(file = delta_pdf, label = "B", title = "Top PBMC vs Ileum pathway differences"),
  list(file = dot_concordant_pdf, label = "C", title = "Top 100 concordant up/down pathways"),
  list(file = heatmap_concordant_pdf, label = "D", title = "Concordant pathway heatmap"),
  list(file = summary_sig_pdf, label = "E", title = "Matched pathways by significance pattern"),
  list(file = summary_concordant_pdf, label = "F", title = "Concordant pathways by cell type")
)

panels <- lapply(panel_info, function(x) {
  make_pdf_panel(x$file, x$label, x$title)
})

title_panel <- ggdraw() +
  draw_label(
    "PBMC versus Ileum immune-cell pathway comparison",
    x = 0.5, y = 0.70,
    hjust = 0.5, vjust = 0.5,
    fontface = "bold", size = 18
  ) +
  draw_label(
    "Integrated multi-panel summary of concordance, difference, and significance across matched Monaco cell types",
    x = 0.5, y = 0.25,
    hjust = 0.5, vjust = 0.5,
    size = 11
  )

top_row <- plot_grid(panels[[1]], panels[[2]], ncol = 2, align = "hv")
mid_row <- plot_grid(panels[[3]], panels[[4]], ncol = 2, align = "hv")
bot_row <- plot_grid(panels[[5]], panels[[6]], ncol = 2, align = "hv")

final_panel <- plot_grid(
  title_panel,
  top_row,
  mid_row,
  bot_row,
  ncol = 1,
  rel_heights = c(0.10, 1, 1, 1)
)

save_plot(final_panel, "PBMC_ILEUM_paper_style_panel.pdf", 16, 22)
save_plot(final_panel, "PBMC_ILEUM_paper_style_panel.png", 16, 22)

cat("Done. Plots written to:", OUT_DIR, "\n")



######