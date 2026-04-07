# =============================================================================
# Cross-Dataset GSEA Exported Significant Pathway Heatmap
#
# Reads:
#   cross_dataset_gsea_output/cross_dataset_GSEA_sig_broad.csv
#   cross_dataset_gsea_output/cross_dataset_GSEA_sig_fine.csv
#
# Makes:
#   - grouped heatmaps
#   - pathway group assignment tables
#
# Grouping:
#   Pathways are grouped into broad biological themes using keyword matching
#   (e.g. IFN/JAK-STAT, Inflammation/TNF, Metabolism, Cell Cycle, etc.)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(fs)
  library(grid)
})

# =============================================================================
# 1. USER SETTINGS
# =============================================================================

OUT_DIR <- "./cross_dataset_gsea_output"
dir_create(OUT_DIR, recurse = TRUE)

SIG_HEATMAP_TOP_N <- 2000
SIG_HEATMAP_USE_ABS_NES <- TRUE
SIG_HEATMAP_SPLIT_BY_DATASET <- TRUE
SIG_HEATMAP_TRUNCATE_ROW_NAMES <- 80

# If TRUE, rows are clustered within each biological theme
SIG_CLUSTER_ROWS <- TRUE

# If TRUE, split columns by dataset
SIG_SPLIT_COLUMNS_BY_DATASET <- TRUE

# =============================================================================
# 2. PATHWAY GROUPING
# =============================================================================

annotate_pathways <- function(pathways) {
  case_when(
    str_detect(pathways, regex("INTERFERON|IFN|JAK|STAT|ISG|IRF", ignore_case = TRUE)) ~ "IFN / JAK-STAT",
    str_detect(pathways, regex("TNF|NFKB|NF.KB|INFLAM|IL1|IL6|CYTOKINE|CHEMOKINE", ignore_case = TRUE)) ~ "Inflammation / TNF",
    str_detect(pathways, regex("ANTIGEN|MHC|HLA|TAP|B2M|CLASS I|CLASS II", ignore_case = TRUE)) ~ "Antigen Presentation",
    str_detect(pathways, regex("T CELL|T-CELL|B CELL|B-CELL|LYMPHOCYTE|LEUKOCYTE|MONOCYTE|MACROPHAGE|NEUTROPHIL", ignore_case = TRUE)) ~ "Immune Cell Programs",
    str_detect(pathways, regex("CELL.?CYCLE|MITOTIC|MITOSIS|E2F|G2M|G2/M|PROLIFERAT", ignore_case = TRUE)) ~ "Cell Cycle / Proliferation",
    str_detect(pathways, regex("APOPTOSIS|DEATH|P53|CASPASE|STRESS", ignore_case = TRUE)) ~ "Cell Death / Stress",
    str_detect(pathways, regex("METABOL|OXIDATIVE|GLYCOL|TCA|RESPIRATORY|MITOCHONDR|FATTY ACID|CHOLESTEROL|LIPID", ignore_case = TRUE)) ~ "Metabolism",
    str_detect(pathways, regex("TRANSLATION|RIBOSOME|PROTEIN SYNTHESIS|MRNA|RNA", ignore_case = TRUE)) ~ "RNA / Translation",
    str_detect(pathways, regex("DNA REPAIR|REPLICATION|CHECKPOINT|CHROMOSOME|DNA DAMAGE", ignore_case = TRUE)) ~ "DNA / Genome Maintenance",
    str_detect(pathways, regex("WNT|MAPK|MTOR|PI3K|AKT|NOTCH|TGF|ERK", ignore_case = TRUE)) ~ "Signaling",
    str_detect(pathways, regex("ADHESION|MIGRATION|CYTOSKELETON|EXTRACELLULAR MATRIX|ECM|INTEGRIN", ignore_case = TRUE)) ~ "Adhesion / Migration",
    TRUE ~ "Other"
  )
}

pathway_group_levels <- c(
  "IFN / JAK-STAT",
  "Inflammation / TNF",
  "Antigen Presentation",
  "Immune Cell Programs",
  "Cell Cycle / Proliferation",
  "Cell Death / Stress",
  "Metabolism",
  "RNA / Translation",
  "DNA / Genome Maintenance",
  "Signaling",
  "Adhesion / Migration",
  "Other"
)

pathway_group_colors <- c(
  "IFN / JAK-STAT" = "#d73027",
  "Inflammation / TNF" = "#fc8d59",
  "Antigen Presentation" = "#4575b4",
  "Immune Cell Programs" = "#74add1",
  "Cell Cycle / Proliferation" = "#984ea3",
  "Cell Death / Stress" = "#f46d43",
  "Metabolism" = "#66bd63",
  "RNA / Translation" = "#1b9e77",
  "DNA / Genome Maintenance" = "#7570b3",
  "Signaling" = "#e6ab02",
  "Adhesion / Migration" = "#a6761d",
  "Other" = "grey70"
)

# =============================================================================
# 3. MAIN FUNCTION
# =============================================================================

plot_exported_sig_heatmap <- function(csv_file,
                                      out_prefix,
                                      top_n = 2000,
                                      use_abs_nes = TRUE,
                                      split_by_dataset = TRUE,
                                      truncate_row_names = 80,
                                      cluster_rows = TRUE) {
  
  if (!file.exists(csv_file)) {
    message(sprintf("[SKIP] File not found: %s", csv_file))
    return(invisible(NULL))
  }
  
  message(sprintf("[READ] %s", csv_file))
  df <- fread(csv_file)
  
  required_cols <- c("dataset", "cluster", "geneset", "Description", "NES", "p.adjust")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing columns in %s: %s",
      csv_file,
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  # ---------------------------------------------------------------------------
  # Clean and define row/column IDs
  # ---------------------------------------------------------------------------
  df <- df %>%
    mutate(
      Description = as.character(Description),
      Description = str_squish(Description),
      row_id = paste0(Description, " [", geneset, "]"),
      col_id = paste(dataset, cluster, sep = " | ")
    )
  
  # Keep strongest NES if duplicate row_id x col_id exists
  df_best <- df %>%
    group_by(row_id, Description, geneset, dataset, cluster, col_id) %>%
    slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # ---------------------------------------------------------------------------
  # Rank pathways globally
  # ---------------------------------------------------------------------------
  pathway_rank <- df_best %>%
    group_by(row_id, Description, geneset) %>%
    summarise(
      max_abs_nes = max(abs(NES), na.rm = TRUE),
      mean_abs_nes = mean(abs(NES), na.rm = TRUE),
      min_fdr = min(p.adjust, na.rm = TRUE),
      n_hits = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(if (use_abs_nes) max_abs_nes else mean_abs_nes), min_fdr)
  
  n_keep <- min(top_n, nrow(pathway_rank))
  
  keep_rows <- pathway_rank %>%
    slice_head(n = n_keep) %>%
    pull(row_id)
  
  df_keep <- df_best %>%
    filter(row_id %in% keep_rows)
  
  if (nrow(df_keep) == 0) {
    message(sprintf("[SKIP] No rows retained for %s", csv_file))
    return(invisible(NULL))
  }
  
  # ---------------------------------------------------------------------------
  # Build matrices
  # ---------------------------------------------------------------------------
  nes_mat <- df_keep %>%
    select(row_id, col_id, NES) %>%
    pivot_wider(names_from = col_id, values_from = NES, values_fill = 0) %>%
    column_to_rownames("row_id") %>%
    as.matrix()
  
  sig_mat <- df_keep %>%
    select(row_id, col_id, p.adjust) %>%
    pivot_wider(names_from = col_id, values_from = p.adjust, values_fill = 1) %>%
    column_to_rownames("row_id") %>%
    as.matrix()
  
  sig_mat <- sig_mat[rownames(nes_mat), colnames(nes_mat), drop = FALSE]
  
  # Preserve ranking order
  row_order_df <- pathway_rank %>%
    filter(row_id %in% rownames(nes_mat)) %>%
    mutate(rank_metric = if (use_abs_nes) max_abs_nes else mean_abs_nes) %>%
    arrange(desc(rank_metric), min_fdr)
  
  nes_mat <- nes_mat[row_order_df$row_id, , drop = FALSE]
  sig_mat <- sig_mat[row_order_df$row_id, , drop = FALSE]
  
  # ---------------------------------------------------------------------------
  # Biological grouping
  # ---------------------------------------------------------------------------
  row_groups <- annotate_pathways(rownames(nes_mat))
  row_groups <- factor(row_groups, levels = pathway_group_levels)
  
  # ---------------------------------------------------------------------------
  # Display row labels
  # ---------------------------------------------------------------------------
  display_rownames <- rownames(nes_mat)
  if (!is.null(truncate_row_names) && truncate_row_names > 0) {
    display_rownames <- str_trunc(
      display_rownames,
      width = truncate_row_names,
      ellipsis = "…"
    )
  }
  
  # ---------------------------------------------------------------------------
  # Column annotation
  # ---------------------------------------------------------------------------
  col_parts <- str_split_fixed(colnames(nes_mat), " \\| ", 2)
  col_dataset <- col_parts[, 1]
  col_cluster <- col_parts[, 2]
  
  dataset_levels <- unique(col_dataset)
  dataset_pal <- setNames(
    colorRampPalette(
      brewer.pal(max(3, min(8, length(dataset_levels))), "Set2")
    )(length(dataset_levels)),
    dataset_levels
  )
  
  top_anno <- HeatmapAnnotation(
    Dataset = col_dataset,
    col = list(Dataset = dataset_pal),
    annotation_name_gp = gpar(fontsize = 9),
    simple_anno_size = unit(4, "mm")
  )
  
  left_anno <- rowAnnotation(
    Category = row_groups,
    col = list(Category = pathway_group_colors),
    annotation_name_gp = gpar(fontsize = 9),
    show_annotation_name = TRUE,
    width = unit(1.2, "cm")
  )
  
  # ---------------------------------------------------------------------------
  # Heatmap colors
  # ---------------------------------------------------------------------------
  max_abs <- max(abs(nes_mat), na.rm = TRUE)
  max_abs <- max(2, min(max_abs, 5))
  
  col_fun <- colorRamp2(
    c(-max_abs, -max_abs / 2, 0, max_abs / 2, max_abs),
    c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D")
  )
  
  column_split <- if (split_by_dataset) {
    factor(col_dataset, levels = unique(col_dataset))
  } else {
    NULL
  }
  
  # ---------------------------------------------------------------------------
  # Build heatmap
  # ---------------------------------------------------------------------------
  ht <- Heatmap(
    nes_mat,
    name = "NES",
    col = col_fun,
    top_annotation = top_anno,
    left_annotation = left_anno,
    row_split = row_groups,
    cluster_rows = cluster_rows,
    cluster_columns = FALSE,
    column_split = column_split,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 8),
    row_labels = display_rownames,
    rect_gp = gpar(col = "white", lwd = 0.3),
    border = TRUE,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    column_title = sprintf(
      "%s exported GSEA heatmap grouped by pathway theme (top %d pathways)",
      out_prefix,
      n_keep
    ),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title = "NES",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(sig_mat[i, j]) && sig_mat[i, j] <= 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 7, fontface = "bold"))
      }
    }
  )
  
  # ---------------------------------------------------------------------------
  # Save outputs
  # ---------------------------------------------------------------------------
  pdf_file <- file.path(
    OUT_DIR,
    paste0(out_prefix, "_exported_sig_heatmap_grouped_top", top_n, ".pdf")
  )
  png_file <- file.path(
    OUT_DIR,
    paste0(out_prefix, "_exported_sig_heatmap_grouped_top", top_n, ".png")
  )
  
  pdf_w <- max(13, ncol(nes_mat) * 0.35 + 5)
  pdf_h <- max(10, nrow(nes_mat) * 0.12 + 4)
  
  pdf(pdf_file, width = pdf_w, height = pdf_h)
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
  dev.off()
  
  png(png_file, width = pdf_w * 150, height = pdf_h * 150, res = 150)
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
  dev.off()
  
  # ---------------------------------------------------------------------------
  # Export group assignment table
  # ---------------------------------------------------------------------------
  grouping_table <- data.frame(
    row_id = rownames(nes_mat),
    pathway_group = as.character(row_groups),
    stringsAsFactors = FALSE
  )
  
  fwrite(
    grouping_table,
    file.path(
      OUT_DIR,
      paste0(out_prefix, "_pathway_group_assignments_top", top_n, ".csv")
    )
  )
  
  message(sprintf("[SAVED] %s", pdf_file))
  message(sprintf("[SAVED] %s", png_file))
  
  invisible(list(
    data = df_keep,
    nes_mat = nes_mat,
    sig_mat = sig_mat,
    row_groups = row_groups,
    heatmap = ht
  ))
}

# =============================================================================
# 4. RUN FOR BROAD AND FINE
# =============================================================================

res_broad_exported <- plot_exported_sig_heatmap(
  csv_file = file.path(OUT_DIR, "cross_dataset_GSEA_sig_broad.csv"),
  out_prefix = "BROAD",
  top_n = SIG_HEATMAP_TOP_N,
  use_abs_nes = SIG_HEATMAP_USE_ABS_NES,
  split_by_dataset = SIG_SPLIT_COLUMNS_BY_DATASET,
  truncate_row_names = SIG_HEATMAP_TRUNCATE_ROW_NAMES,
  cluster_rows = SIG_CLUSTER_ROWS
)

res_fine_exported <- plot_exported_sig_heatmap(
  csv_file = file.path(OUT_DIR, "cross_dataset_GSEA_sig_fine.csv"),
  out_prefix = "FINE",
  top_n = SIG_HEATMAP_TOP_N,
  use_abs_nes = SIG_HEATMAP_USE_ABS_NES,
  split_by_dataset = SIG_SPLIT_COLUMNS_BY_DATASET,
  truncate_row_names = SIG_HEATMAP_TRUNCATE_ROW_NAMES,
  cluster_rows = SIG_CLUSTER_ROWS
)

# =============================================================================
# 5. OPTIONAL COMBINED TABLE OF GROUP COUNTS
# =============================================================================

combine_group_counts <- function(res_obj, label) {
  if (is.null(res_obj)) return(NULL)
  data.frame(
    analysis = label,
    pathway_group = as.character(res_obj$row_groups),
    stringsAsFactors = FALSE
  ) %>%
    count(analysis, pathway_group, name = "n_pathways")
}

group_count_table <- bind_rows(
  combine_group_counts(res_broad_exported, "BROAD"),
  combine_group_counts(res_fine_exported, "FINE")
)

if (nrow(group_count_table) > 0) {
  fwrite(
    group_count_table,
    file.path(OUT_DIR, "grouped_pathway_counts_summary.csv")
  )
}

message("=== Done ===")
message(sprintf("Outputs written to: %s", OUT_DIR))