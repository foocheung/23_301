# =========================
# TOP CONCORDANT + DISCORDANT HEATMAPS
# PBMC and DUODENUM paired together by cell type
# =========================

# -------------------------
# Helper: order columns as PBMC then DUODENUM within each cell type
# -------------------------
pair_columns_by_celltype <- function(mat) {
  cell_types <- sub("^(PBMC|DUODENUM) \\| ", "", colnames(mat))
  cell_types_unique <- unique(cell_types)
  
  desired_col_order <- unlist(lapply(cell_types_unique, function(ct) {
    c(
      paste0("PBMC | ", ct),
      paste0("DUODENUM | ", ct)
    )
  }))
  
  desired_col_order <- desired_col_order[desired_col_order %in% colnames(mat)]
  mat <- mat[, desired_col_order, drop = FALSE]
  
  column_split_vec <- sub("^(PBMC|DUODENUM) \\| ", "", colnames(mat))
  column_split_vec <- factor(column_split_vec, levels = unique(column_split_vec))
  
  list(
    mat = mat,
    column_split_vec = column_split_vec
  )
}

# -------------------------
# CONCORDANT
# Top Both Up + Top Both Down
# -------------------------
top_concordant <- matched %>%
  filter(
    sig_status != "Neither significant",
    quadrant %in% c(
      "Q1: Both Up",
      "Q3: Both Down"
    )
  ) %>%
  mutate(
    direction_group = case_when(
      quadrant == "Q1: Both Up"   ~ "Both Up",
      quadrant == "Q3: Both Down" ~ "Both Down"
    ),
    concordance_score = abs(NES__PBMC) + abs(NES__ILEUM)
  ) %>%
  group_by(direction_group) %>%
  slice_max(order_by = concordance_score, n = 50, with_ties = FALSE) %>%
  ungroup()

concordant_heatmap_long <- make_dual_dataset_long(top_concordant, include_direction = TRUE)
concordant_heatmap_mat  <- build_nes_matrix(concordant_heatmap_long)

row_split_df_conc <- top_concordant %>%
  distinct(pathway_id, Description, geneset, direction_group) %>%
  mutate(row_id = paste0(Description, " [", geneset, "]")) %>%
  distinct(row_id, direction_group)

row_split_vec_conc <- row_split_df_conc$direction_group[
  match(rownames(concordant_heatmap_mat), row_split_df_conc$row_id)
]

row_split_vec_conc <- factor(
  row_split_vec_conc,
  levels = c("Both Up", "Both Down")
)

conc_cols <- pair_columns_by_celltype(concordant_heatmap_mat)
concordant_heatmap_mat <- conc_cols$mat
col_split_vec_conc <- conc_cols$column_split_vec

heatmap_concordant_pdf <- save_heatmap_pdf(
  mat             = concordant_heatmap_mat,
  filename        = "heatmap_top10_both_up_and_down.pdf",
  title           = "Top pathways Both Up and Top pathways Both Down\n(PBMC and Duodenum paired within each cell type)",
  row_split       = row_split_vec_conc,
  column_split    = col_split_vec_conc,
  cluster_columns = FALSE
)

# -------------------------
# DISCORDANT
# Top PBMC Up / DUODENUM Down + Top DUODENUM Up / PBMC Down
# -------------------------
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
    discordance_score = abs(NES__PBMC) + abs(NES__ILEUM)
  ) %>%
  group_by(direction_group) %>%
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

disc_cols <- pair_columns_by_celltype(discordant_heatmap_mat)
discordant_heatmap_mat <- disc_cols$mat
col_split_vec_disc <- disc_cols$column_split_vec

heatmap_discordant_pdf <- save_heatmap_pdf(
  mat             = discordant_heatmap_mat,
  filename        = "heatmap_top10_discordant_paired.pdf",
  title           = "Top pathways PBMC Up / Duodenum Down and\n Top pathways Duodenum Up / PBMC Down\n(PBMC and Duodenum paired within each cell type)",
  row_split       = row_split_vec_disc,
  column_split    = col_split_vec_disc,
  cluster_columns = FALSE
)

# -------------------------
# Quick checks
# -------------------------
cat("Concordant rows:", nrow(top_concordant), "\n")
cat("Discordant rows:", nrow(top_discordant), "\n")

cat("Concordant columns:\n")
print(colnames(concordant_heatmap_mat))
print(col_split_vec_conc)

cat("Discordant columns:\n")
print(colnames(discordant_heatmap_mat))
print(col_split_vec_disc)