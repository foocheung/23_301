###############################################################################
# PLOT: MildTreated vs MildUntreated (paired) — across 3 annotation strategies
#
# Reads DEG_side_by_side_summary.csv and GSEA_side_by_side_summary.csv from:
#   - deg_gsea_pooled_vs_within_pool_Monaco_main_pruned/
#   - deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned/
#   - deg_gsea_pooled_vs_within_pool_cluster/
#
# Produces:
#   1)  DEG_side_by_side_total.png
#   2)  DEG_side_by_side_up_down.png
#   3)  GSEA_side_by_side_by_database.png
#   4)  GSEA_side_by_side_total.png
#   5)  DEG_heatmap_annotation_overlap.png   (n_sig per celltype x annotation)
#   6)  GSEA_heatmap_annotation_overlap.png  (n_sig pathways per celltype x annotation)
#   7)  DEG_dotplot_up_down_ratio.png        (logFC direction bias by annotation)
#   8)  GSEA_upset_shared_pathways.png       (shared significant pathways across annotations)
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

# =============================================================================
# SETTINGS — edit these paths if needed
# =============================================================================
BASE_DIRS <- list(
  Monaco_main  = "deg_gsea_pooled_vs_within_pool_Monaco_main_pruned",
  Monaco_fine  = "deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned",
  Cluster      = "deg_gsea_pooled_vs_within_pool_cluster"
)

OUT_DIR <- "paired_comparison_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAIRED_CONTRAST <- "MildTreated_vs_MildUntreated_paired"

# =============================================================================
# HELPERS
# =============================================================================



save_plot <- function(p, path_no_ext, width = 14, height = 8) {
  ggplot2::ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
  message("Saved: ", path_no_ext)
}

annotation_palette <- c(
  Monaco_main = "#2166AC",
  Monaco_fine = "#D6604D",
  Cluster     = "#4DAC26"
)

annotation_shapes <- c(
  Monaco_main = 16,
  Monaco_fine = 17,
  Cluster     = 15
)

# =============================================================================
# LOAD DATA
# =============================================================================
load_summaries <- function(base_dirs, contrast_filter) {
  deg_list  <- list()
  gsea_list <- list()
  
  standardize_df <- function(df, ann_name) {
    if ("celltype" %in% names(df)) {
      df$celltype <- as.character(df$celltype)
    }
    df$annotation <- ann_name
    df
  }
  
  for (ann_name in names(base_dirs)) {
    bdir <- base_dirs[[ann_name]]
    
    deg_path  <- file.path(bdir, "DEG_side_by_side_summary.csv")
    gsea_path <- file.path(bdir, "GSEA_side_by_side_summary.csv")
    
    if (file.exists(deg_path)) {
      deg_list[[ann_name]] <- data.table::fread(deg_path) %>%
        dplyr::filter(contrast == contrast_filter) %>%
        standardize_df(ann_name)
    } else {
      warning("DEG summary not found: ", deg_path)
    }
    
    if (file.exists(gsea_path)) {
      gsea_list[[ann_name]] <- data.table::fread(gsea_path) %>%
        dplyr::filter(contrast == contrast_filter) %>%
        standardize_df(ann_name)
    } else {
      warning("GSEA summary not found: ", gsea_path)
    }
  }
  
  list(
    deg  = dplyr::bind_rows(deg_list),
    gsea = dplyr::bind_rows(gsea_list)
  )
}

message("Loading paired contrast data...")
dat <- load_summaries(BASE_DIRS, PAIRED_CONTRAST)
deg_df  <- dat$deg
gsea_df <- dat$gsea

if (nrow(deg_df) == 0 && nrow(gsea_df) == 0) {
  stop(
    "No data found for contrast '", PAIRED_CONTRAST, "'. ",
    "Check that:\n",
    "  1) The base directory paths are correct.\n",
    "  2) The paired analysis ran successfully.\n",
    "  3) DEG_side_by_side_summary.csv exists in each output directory."
  )
}

message("DEG rows loaded:  ", nrow(deg_df))
message("GSEA rows loaded: ", nrow(gsea_df))

# =============================================================================
# THEME
# =============================================================================
theme_paired <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "#F0F0F0", color = NA),
      strip.text        = ggplot2::element_text(face = "bold", size = base_size - 1),
      axis.text.x       = ggplot2::element_text(angle = 35, hjust = 1, size = base_size - 2),
      axis.text.y       = ggplot2::element_text(size = base_size - 2),
      legend.position   = "right",
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle     = ggplot2::element_text(color = "grey40", size = base_size - 1)
    )
}

contrast_label <- "Mild Treated vs Mild Untreated (paired)"

# =============================================================================
# PLOT 1: DEG_side_by_side_total
# n_sig per celltype, grouped by annotation
# =============================================================================
if (nrow(deg_df) > 0) {
  message("Plotting DEG_side_by_side_total...")
  
  p1 <- ggplot2::ggplot(
    deg_df,
    ggplot2::aes(
      x    = forcats::fct_reorder(celltype, n_sig, .fun = max, .desc = TRUE),
      y    = n_sig,
      fill = annotation
    )
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_sig),
      position = ggplot2::position_dodge(width = 0.8),
      vjust = -0.35,
      size  = 3
    ) +
    ggplot2::scale_fill_manual(
      values = annotation_palette,
      name   = "Annotation\nstrategy"
    ) +
    ggplot2::labs(
      title    = "Total significant DEGs — Paired contrast",
      subtitle = contrast_label,
      x        = "Cell type",
      y        = "Number of significant DEGs\n(FDR < 0.05, |logFC| ≥ 0.25)"
    ) +
    theme_paired()
  
  save_plot(p1, file.path(OUT_DIR, "DEG_side_by_side_total"), width = 16, height = 7)
}

# =============================================================================
# PLOT 2: DEG_side_by_side_up_down
# Stacked up/down bars per celltype, faceted by annotation
# =============================================================================
if (nrow(deg_df) > 0) {
  message("Plotting DEG_side_by_side_up_down...")
  
  deg_long <- deg_df %>%
    tidyr::pivot_longer(
      cols      = c(n_up, n_down),
      names_to  = "direction",
      values_to = "n_genes"
    ) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        direction == "n_up"   ~ "Up-regulated",
        direction == "n_down" ~ "Down-regulated"
      ),
      direction = factor(direction, levels = c("Up-regulated", "Down-regulated")),
      celltype  = forcats::fct_reorder(celltype, n_sig, .fun = max, .desc = TRUE)
    )
  
  p2 <- ggplot2::ggplot(
    deg_long,
    ggplot2::aes(x = celltype, y = n_genes, fill = direction)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(n_genes > 0, n_genes, "")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 2.8,
      color = "white",
      fontface = "bold"
    ) +
    ggplot2::facet_wrap(~ annotation, ncol = 3) +
    ggplot2::scale_fill_manual(
      values = c("Up-regulated" = "#C0392B", "Down-regulated" = "#2980B9"),
      name   = "Direction"
    ) +
    ggplot2::labs(
      title    = "Up- and down-regulated DEGs — Paired contrast",
      subtitle = contrast_label,
      x        = "Cell type",
      y        = "Number of genes"
    ) +
    theme_paired()
  
  save_plot(p2, file.path(OUT_DIR, "DEG_side_by_side_up_down"), width = 18, height = 7)
}

# =============================================================================
# PLOT 3: GSEA_side_by_side_by_database
# Significant pathways per database, faceted by celltype × annotation
# =============================================================================
if (nrow(gsea_df) > 0) {
  message("Plotting GSEA_side_by_side_by_database...")
  
  gsea_ann <- gsea_df %>%
    dplyr::mutate(
      celltype = forcats::fct_reorder(celltype, n_sig, .fun = sum, .desc = TRUE)
    )
  
  p3 <- ggplot2::ggplot(
    gsea_ann,
    ggplot2::aes(x = database, y = n_sig, fill = database)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_sig),
      vjust = -0.3,
      size  = 2.8
    ) +
    ggplot2::facet_grid(annotation ~ celltype, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c(
        GO_BP    = "#1B7837",
        Reactome = "#762A83",
        Hallmark = "#E08214"
      ),
      name = "Database"
    ) +
    ggplot2::labs(
      title    = "Significant GSEA pathways by database — Paired contrast",
      subtitle = contrast_label,
      x        = "Database",
      y        = "Number of significant pathways (padj < 0.05)"
    ) +
    theme_paired() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "bottom"
    )
  
  save_plot(p3, file.path(OUT_DIR, "GSEA_side_by_side_by_database"), width = 22, height = 10)
}

# =============================================================================
# PLOT 4: GSEA_side_by_side_total
# Total significant pathways per celltype, stacked by database, faceted by annotation
# =============================================================================
if (nrow(gsea_df) > 0) {
  message("Plotting GSEA_side_by_side_total...")
  
  gsea_total <- gsea_df %>%
    dplyr::group_by(celltype, annotation) %>%
    dplyr::summarise(n_sig = sum(n_sig, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      celltype = forcats::fct_reorder(celltype, n_sig, .fun = max, .desc = TRUE)
    )
  
  p4 <- ggplot2::ggplot(
    gsea_total,
    ggplot2::aes(
      x    = celltype,
      y    = n_sig,
      fill = annotation
    )
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_sig),
      position = ggplot2::position_dodge(width = 0.8),
      vjust = -0.35,
      size  = 3
    ) +
    ggplot2::scale_fill_manual(
      values = annotation_palette,
      name   = "Annotation\nstrategy"
    ) +
    ggplot2::labs(
      title    = "Total significant GSEA pathways — Paired contrast",
      subtitle = paste0(contrast_label, " | all databases combined"),
      x        = "Cell type",
      y        = "Total significant pathways (padj < 0.05)"
    ) +
    theme_paired()
  
  save_plot(p4, file.path(OUT_DIR, "GSEA_side_by_side_total"), width = 16, height = 7)
}

# =============================================================================
# PLOT 5 (EXTRA): Heatmap — DEG n_sig per celltype × annotation
# =============================================================================
if (nrow(deg_df) > 0) {
  message("Plotting DEG_heatmap_annotation_overlap...")
  
  deg_heat <- deg_df %>%
    dplyr::select(celltype, annotation, n_sig) %>%
    dplyr::mutate(
      celltype   = forcats::fct_reorder(celltype, n_sig, .fun = sum, .desc = TRUE),
      annotation = factor(annotation, levels = names(BASE_DIRS))
    )
  
  p5 <- ggplot2::ggplot(
    deg_heat,
    ggplot2::aes(x = annotation, y = celltype, fill = n_sig)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_sig, color = n_sig > max(n_sig) * 0.6),
      size = 3.5, fontface = "bold"
    ) +
    ggplot2::scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
    ggplot2::scale_fill_gradientn(
      colors = c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"),
      name   = "Sig. DEGs",
      na.value = "grey90"
    ) +
    ggplot2::labs(
      title    = "DEG heatmap: significant genes per cell type & annotation",
      subtitle = contrast_label,
      x        = "Annotation strategy",
      y        = "Cell type"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      plot.title  = ggplot2::element_text(face = "bold")
    )
  
  save_plot(p5, file.path(OUT_DIR, "DEG_heatmap_annotation_overlap"), width = 10, height = 9)
}

# =============================================================================
# PLOT 6 (EXTRA): Heatmap — GSEA n_sig per celltype × annotation (per database)
# =============================================================================
if (nrow(gsea_df) > 0) {
  message("Plotting GSEA_heatmap_annotation_overlap...")
  
  gsea_heat <- gsea_df %>%
    dplyr::mutate(
      annotation = factor(annotation, levels = names(BASE_DIRS)),
      celltype   = forcats::fct_reorder(celltype, n_sig, .fun = sum, .desc = TRUE)
    )
  
  p6 <- ggplot2::ggplot(
    gsea_heat,
    ggplot2::aes(x = annotation, y = celltype, fill = n_sig)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_sig, color = n_sig > max(n_sig) * 0.6),
      size = 3
    ) +
    ggplot2::scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
    ggplot2::scale_fill_gradientn(
      colors   = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FC8D59", "#B30000"),
      name     = "Sig. pathways",
      na.value = "grey90"
    ) +
    ggplot2::facet_wrap(~ database, ncol = 3) +
    ggplot2::labs(
      title    = "GSEA heatmap: significant pathways per cell type, annotation & database",
      subtitle = contrast_label,
      x        = "Annotation strategy",
      y        = "Cell type"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#F0F0F0"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 20, hjust = 1),
      plot.title       = ggplot2::element_text(face = "bold")
    )
  
  save_plot(p6, file.path(OUT_DIR, "GSEA_heatmap_annotation_overlap"), width = 18, height = 10)
}

# =============================================================================
# PLOT 7 (EXTRA): DEG Direction bias — up/down ratio dotplot per annotation
# =============================================================================
if (nrow(deg_df) > 0) {
  message("Plotting DEG_dotplot_up_down_ratio...")
  
  # deg_ratio <- deg_df %>%
  #   dplyr::mutate(
  #     total       = n_up + n_down,
  #     pct_up      = ifelse(total > 0, 100 * n_up / total, NA_real_),
  #     annotation  = factor(annotation, levels = names(BASE_DIRS)),
  #     celltype    = forcats::fct_reorder(celltype, pct_up, .fun = mean, na.rm = TRUE, .desc = TRUE)
  #   ) %>%
  #   dplyr::filter(!is.na(pct_up), total > 0)
  # 
 
    deg_ratio <- deg_df %>%
    dplyr::mutate(
      total      = n_up + n_down,
      pct_up     = dplyr::if_else(total > 0, 100 * n_up / total, NA_real_),
      annotation = factor(annotation, levels = names(BASE_DIRS)),
      celltype   = as.character(celltype)
    ) %>%
    dplyr::filter(total > 0, !is.na(pct_up), is.finite(pct_up))
  
  if (nrow(deg_ratio) > 0) {
    celltype_order <- deg_ratio %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(mean_pct_up = mean(pct_up, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(mean_pct_up)) %>%
      dplyr::pull(celltype)
    
    deg_ratio <- deg_ratio %>%
      dplyr::mutate(
        celltype = factor(celltype, levels = celltype_order)
      )
  }
  
  if (nrow(deg_ratio) > 0) {
    p7 <- ggplot2::ggplot(
      deg_ratio,
      ggplot2::aes(
        x     = pct_up,
        y     = celltype,
        color = annotation,
        size  = total,
        shape = annotation
      )
    ) +
      ggplot2::geom_vline(xintercept = 50, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_color_manual(
        values = annotation_palette,
        name   = "Annotation"
      ) +
      ggplot2::scale_shape_manual(
        values = annotation_shapes,
        name   = "Annotation"
      ) +
      ggplot2::scale_size_continuous(
        name   = "Total sig. DEGs",
        range  = c(2, 10),
        breaks = c(5, 20, 50, 100, 200)
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 100),
        labels = scales::label_percent(scale = 1),
        name   = "% up-regulated (of significant DEGs)"
      ) +
      ggplot2::labs(
        title    = "Direction bias: % up-regulated DEGs per cell type & annotation",
        subtitle = paste0(contrast_label, "\nDot size = total significant DEGs; dashed line = 50%"),
        y        = "Cell type"
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      )
    
    save_plot(p7, file.path(OUT_DIR, "DEG_dotplot_up_down_ratio"), width = 12, height = 9)
  } else {
    message("Skipping DEG_dotplot_up_down_ratio: no rows with n_up + n_down > 0")
  }
}

# =============================================================================
# PLOT 8 (EXTRA): Concordance scatter — Monaco_main vs Cluster n_sig DEGs
# For shared cell types across annotation strategies
# =============================================================================
if (nrow(deg_df) > 0 && length(unique(deg_df$annotation)) >= 2) {
  message("Plotting DEG_concordance_scatter...")
  
  ann_pairs <- combn(names(BASE_DIRS), 2, simplify = FALSE)
  
  concordance_plots <- lapply(ann_pairs, function(pair) {
    a1 <- pair[1]; a2 <- pair[2]
    
    d1 <- deg_df %>% dplyr::filter(annotation == a1) %>%
      dplyr::select(celltype, n_sig_a1 = n_sig, n_up_a1 = n_up, n_down_a1 = n_down)
    d2 <- deg_df %>% dplyr::filter(annotation == a2) %>%
      dplyr::select(celltype, n_sig_a2 = n_sig, n_up_a2 = n_up, n_down_a2 = n_down)
    
    both <- dplyr::inner_join(d1, d2, by = "celltype")
    if (nrow(both) < 2) return(NULL)
    
    cor_val <- round(cor(both$n_sig_a1, both$n_sig_a2, method = "spearman", use = "complete.obs"), 3)
    
    ggplot2::ggplot(both, ggplot2::aes(x = n_sig_a1, y = n_sig_a2, label = celltype)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "dashed") +
      ggplot2::geom_point(size = 3, color = annotation_palette[a1], alpha = 0.8) +
      ggplot2::geom_text(
        vjust = -0.6, size = 2.8, color = "grey30",
        check_overlap = TRUE
      ) +
      ggplot2::annotate(
        "text", x = -Inf, y = Inf,
        label = paste0("ρ = ", cor_val),
        hjust = -0.1, vjust = 1.5, size = 4.5, color = "grey20", fontface = "bold"
      ) +
      ggplot2::labs(
        x = paste0("Sig. DEGs (", a1, ")"),
        y = paste0("Sig. DEGs (", a2, ")"),
        title = paste0(a1, " vs ", a2)
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  })
  
  concordance_plots <- Filter(Negate(is.null), concordance_plots)
  
  if (length(concordance_plots) > 0) {
    p8 <- patchwork::wrap_plots(concordance_plots, ncol = length(concordance_plots)) +
      patchwork::plot_annotation(
        title    = "DEG concordance across annotation strategies — Paired contrast",
        subtitle = paste0(contrast_label, "\nSpearman ρ computed on n_sig per shared cell type"),
        theme    = ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(color = "grey40", size = 11)
        )
      )
    
    save_plot(p8, file.path(OUT_DIR, "DEG_concordance_scatter"), width = 16, height = 6)
  }
}

# =============================================================================
# PLOT 9 (EXTRA): GSEA concordance scatter per database
# =============================================================================
if (nrow(gsea_df) > 0 && length(unique(gsea_df$annotation)) >= 2) {
  message("Plotting GSEA_concordance_scatter...")
  
  ann_pairs <- combn(names(BASE_DIRS), 2, simplify = FALSE)
  databases <- unique(gsea_df$database)
  
  conc_list <- list()
  
  for (db in databases) {
    for (pair in ann_pairs) {
      a1 <- pair[1]; a2 <- pair[2]
      
      d1 <- gsea_df %>%
        dplyr::filter(annotation == a1, database == db) %>%
        dplyr::select(celltype, n_sig_a1 = n_sig)
      d2 <- gsea_df %>%
        dplyr::filter(annotation == a2, database == db) %>%
        dplyr::select(celltype, n_sig_a2 = n_sig)
      
      both <- dplyr::inner_join(d1, d2, by = "celltype")
      if (nrow(both) < 2) next
      
      cor_val <- round(cor(both$n_sig_a1, both$n_sig_a2, method = "spearman", use = "complete.obs"), 3)
      
      g <- ggplot2::ggplot(both, ggplot2::aes(x = n_sig_a1, y = n_sig_a2, label = celltype)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "dashed") +
        ggplot2::geom_point(size = 3, color = annotation_palette[a1], alpha = 0.8) +
        ggplot2::geom_text(vjust = -0.6, size = 2.5, check_overlap = TRUE) +
        ggplot2::annotate(
          "text", x = -Inf, y = Inf,
          label = paste0("ρ = ", cor_val),
          hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold"
        ) +
        ggplot2::labs(
          x     = paste0(a1),
          y     = paste0(a2),
          title = paste0(db, ": ", a1, " vs ", a2)
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 9))
      
      conc_list[[length(conc_list) + 1]] <- g
    }
  }
  
  if (length(conc_list) > 0) {
    ncols <- length(ann_pairs)
    p9 <- patchwork::wrap_plots(conc_list, ncol = ncols) +
      patchwork::plot_annotation(
        title    = "GSEA pathway concordance across annotation strategies — Paired contrast",
        subtitle = paste0(contrast_label, "\nRows = databases, columns = annotation pairs"),
        theme    = ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(color = "grey40", size = 11)
        )
      )
    
    save_plot(p9, file.path(OUT_DIR, "GSEA_concordance_scatter"), width = 18, height = 6 * ceiling(length(databases)))
  }
}


# =============================================================================
# EXPORT: spot-check workbook of top positive/negative DEGs and pathways
# =============================================================================
message("\nExporting spot-check workbook...")
spotcheck <- export_spotcheck_workbook(
  base_dirs = BASE_DIRS,
  contrast_filter = PAIRED_CONTRAST,
  out_dir = OUT_DIR,
  top_n = 10
)


# =============================================================================
# SUMMARY TABLE
# =============================================================================
message("\n=== Summary: paired contrast across annotation strategies ===")

if (nrow(deg_df) > 0) {
  cat("\n--- DEG summary ---\n")
  print(
    deg_df %>%
      dplyr::select(annotation, celltype, n_up, n_down, n_sig) %>%
      dplyr::arrange(annotation, dplyr::desc(n_sig))
  )
}

if (nrow(gsea_df) > 0) {
  cat("\n--- GSEA summary ---\n")
  print(
    gsea_df %>%
      dplyr::group_by(annotation, celltype) %>%
      dplyr::summarise(total_pathways = sum(n_sig), .groups = "drop") %>%
      dplyr::arrange(annotation, dplyr::desc(total_pathways))
  )
}

message("\nAll plots saved to: ", OUT_DIR)
message("Done.")


#######
##############
#####################################3

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

# =============================================================================
# SETTINGS
# =============================================================================
BASE_DIRS <- list(
  Monaco_main  = "deg_gsea_pooled_vs_within_pool_Monaco_main_pruned",
  Monaco_fine  = "deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned",
  Cluster      = "deg_gsea_pooled_vs_within_pool_cluster"
)

OUT_DIR <- "paired_comparison_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAIRED_CONTRAST <- "MildTreated_vs_MildUntreated_paired"

# =============================================================================
# HELPERS
# =============================================================================
first_existing_col <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

safe_read_csv <- function(path) {
  tryCatch(
    data.table::fread(path),
    error = function(e) NULL
  )
}

extract_celltype_from_path <- function(path) {
  parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/")[[1]]
  parts_rev <- rev(parts)
  
  ignore_parts <- c(
    "deg_results", "gsea_results", "DEG", "GSEA",
    "GO", "GO_BP", "Hallmark", "Reactome",
    "paired_cross_pool", "pooled", "gsea", "deg"
  )
  
  for (p in parts_rev) {
    if (!grepl("\\.csv$", p, ignore.case = TRUE) && !p %in% ignore_parts) {
      return(p)
    }
  }
  
  return(NA_character_)
}

looks_like_deg_file <- function(df, file_path) {
  has_gene <- any(c("gene", "Gene", "symbol", "SYMBOL", "feature", "rowname") %in% colnames(df))
  has_logfc <- any(c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC") %in% colnames(df))
  has_padj <- any(c("padj", "FDR", "adj.P.Val", "p_val_adj") %in% colnames(df))
  name_hint <- grepl("deg|marker|deseq|limma", basename(file_path), ignore.case = TRUE)
  (has_gene && has_logfc && has_padj) || (name_hint && has_gene && has_padj)
}

looks_like_gsea_file <- function(df, file_path) {
  has_pathway <- any(c("Description", "pathway", "Pathway", "term", "ID") %in% colnames(df))
  has_nes <- any(c("NES", "nes") %in% colnames(df))
  has_padj <- any(c("padj", "p.adjust", "FDR", "qvalue") %in% colnames(df))
  name_hint <- grepl("gsea|hallmark|reactome|go", basename(file_path), ignore.case = TRUE)
  (has_pathway && has_nes && has_padj) || (name_hint && has_pathway && has_padj)
}

file_matches_contrast <- function(df, file_path, contrast_filter) {
  contrast_col <- first_existing_col(df, c("contrast", "comparison", "Contrast"))
  
  if (!is.na(contrast_col)) {
    vals <- unique(as.character(df[[contrast_col]]))
    return(any(vals == contrast_filter))
  }
  
  file_text <- paste(dirname(file_path), basename(file_path), sep = "/")
  
  exact_hit <- grepl(contrast_filter, file_text, fixed = TRUE)
  loose_hit <- grepl("MildTreated", file_text, ignore.case = TRUE) &&
    grepl("MildUntreated", file_text, ignore.case = TRUE)
  
  exact_hit || loose_hit
}

inspect_result_files <- function(base_dirs, out_dir) {
  out_list <- list()
  
  for (ann_name in names(base_dirs)) {
    bdir <- base_dirs[[ann_name]]
    
    csvs <- list.files(
      bdir,
      pattern = "\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(csvs) == 0) next
    
    ann_info <- lapply(csvs, function(f) {
      df <- safe_read_csv(f)
      if (is.null(df)) {
        return(data.frame(
          annotation = ann_name,
          file = f,
          file_name = basename(f),
          detected_type = "Unreadable",
          n_rows = NA_integer_,
          columns = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
      
      detected_type <- dplyr::case_when(
        looks_like_deg_file(df, f)  ~ "DEG_like",
        looks_like_gsea_file(df, f) ~ "GSEA_like",
        TRUE                        ~ "Other"
      )
      
      data.frame(
        annotation = ann_name,
        file = f,
        file_name = basename(f),
        detected_type = detected_type,
        n_rows = nrow(df),
        columns = paste(colnames(df), collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
    
    out_list[[ann_name]] <- dplyr::bind_rows(ann_info)
  }
  
  inventory <- dplyr::bind_rows(out_list)
  
  if (nrow(inventory) > 0) {
    data.table::fwrite(
      inventory,
      file.path(out_dir, "result_file_inventory.csv")
    )
    message("Saved file inventory: ", file.path(out_dir, "result_file_inventory.csv"))
  }
  
  invisible(inventory)
}

collect_detailed_deg_tables <- function(base_dirs, contrast_filter, top_n = 10) {
  out_list <- list()
  
  for (ann_name in names(base_dirs)) {
    bdir <- base_dirs[[ann_name]]
    
    csvs <- list.files(
      bdir,
      pattern = "\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    csvs <- csvs[!basename(csvs) %in% c("DEG_side_by_side_summary.csv", "GSEA_side_by_side_summary.csv")]
    
    for (f in csvs) {
      df <- safe_read_csv(f)
      if (is.null(df) || nrow(df) == 0) next
      if (!looks_like_deg_file(df, f)) next
      if (!file_matches_contrast(df, f, contrast_filter)) next
      
      gene_col  <- first_existing_col(df, c("gene", "Gene", "symbol", "SYMBOL", "feature", "rowname"))
      logfc_col <- first_existing_col(df, c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC"))
      padj_col  <- first_existing_col(df, c("padj", "FDR", "adj.P.Val", "p_val_adj"))
      cell_col  <- first_existing_col(df, c("celltype", "cell_type", "cluster", "Cluster"))
      
      if (is.na(gene_col) || is.na(logfc_col) || is.na(padj_col)) next
      
      tmp <- df
      
      contrast_col <- first_existing_col(tmp, c("contrast", "comparison", "Contrast"))
      if (!is.na(contrast_col)) {
        tmp <- tmp %>% dplyr::filter(as.character(.data[[contrast_col]]) == contrast_filter)
      }
      if (nrow(tmp) == 0) next
      
      tmp <- tmp %>%
        dplyr::mutate(
          annotation  = ann_name,
          source_file = basename(f),
          source_path = f,
          celltype = if (!is.na(cell_col)) as.character(.data[[cell_col]]) else extract_celltype_from_path(f),
          gene     = as.character(.data[[gene_col]]),
          logFC    = suppressWarnings(as.numeric(.data[[logfc_col]])),
          padj     = suppressWarnings(as.numeric(.data[[padj_col]]))
        ) %>%
        dplyr::filter(
          !is.na(celltype), celltype != "",
          !is.na(gene), gene != "",
          !is.na(logFC),
          !is.na(padj),
          is.finite(logFC),
          is.finite(padj),
          padj < 0.05
        )
      
      if (nrow(tmp) == 0) next
      
      pos <- tmp %>%
        dplyr::filter(logFC > 0) %>%
        dplyr::group_by(annotation, celltype) %>%
        dplyr::arrange(dplyr::desc(logFC), padj, .by_group = TRUE) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Positive")
      
      neg <- tmp %>%
        dplyr::filter(logFC < 0) %>%
        dplyr::group_by(annotation, celltype) %>%
        dplyr::arrange(logFC, padj, .by_group = TRUE) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Negative")
      
      picked <- dplyr::bind_rows(pos, neg) %>%
        dplyr::ungroup() %>%
        dplyr::select(annotation, celltype, direction, gene, logFC, padj, source_file, source_path)
      
      if (nrow(picked) > 0) {
        out_list[[paste0(ann_name, "::", f)]] <- picked
      }
    }
  }
  
  if (length(out_list) == 0) return(data.frame())
  
  dplyr::bind_rows(out_list) %>%
    dplyr::group_by(annotation, celltype, direction, gene) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(annotation, celltype, direction) %>%
    dplyr::arrange(
      dplyr::if_else(direction == "Positive", -logFC, logFC),
      padj,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
}

collect_detailed_gsea_tables <- function(base_dirs, contrast_filter, top_n = 10) {
  out_list <- list()
  
  for (ann_name in names(base_dirs)) {
    bdir <- base_dirs[[ann_name]]
    
    csvs <- list.files(
      bdir,
      pattern = "\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    csvs <- csvs[!basename(csvs) %in% c("DEG_side_by_side_summary.csv", "GSEA_side_by_side_summary.csv")]
    
    for (f in csvs) {
      df <- safe_read_csv(f)
      if (is.null(df) || nrow(df) == 0) next
      if (!looks_like_gsea_file(df, f)) next
      if (!file_matches_contrast(df, f, contrast_filter)) next
      
      pathway_col <- first_existing_col(df, c("Description", "pathway", "Pathway", "term", "ID"))
      nes_col     <- first_existing_col(df, c("NES", "nes"))
      padj_col    <- first_existing_col(df, c("padj", "p.adjust", "FDR", "qvalue"))
      cell_col    <- first_existing_col(df, c("celltype", "cell_type", "cluster", "Cluster"))
      db_col      <- first_existing_col(df, c("database", "geneset", "collection"))
      
      if (is.na(pathway_col) || is.na(nes_col) || is.na(padj_col)) next
      
      tmp <- df
      
      contrast_col <- first_existing_col(tmp, c("contrast", "comparison", "Contrast"))
      if (!is.na(contrast_col)) {
        tmp <- tmp %>% dplyr::filter(as.character(.data[[contrast_col]]) == contrast_filter)
      }
      if (nrow(tmp) == 0) next
      
      tmp <- tmp %>%
        dplyr::mutate(
          annotation  = ann_name,
          source_file = basename(f),
          source_path = f,
          celltype = if (!is.na(cell_col)) as.character(.data[[cell_col]]) else extract_celltype_from_path(f),
          pathway  = as.character(.data[[pathway_col]]),
          NES      = suppressWarnings(as.numeric(.data[[nes_col]])),
          padj     = suppressWarnings(as.numeric(.data[[padj_col]])),
          database = if (!is.na(db_col)) as.character(.data[[db_col]]) else "Unknown"
        ) %>%
        dplyr::filter(
          !is.na(celltype), celltype != "",
          !is.na(pathway), pathway != "",
          !is.na(NES),
          !is.na(padj),
          is.finite(NES),
          is.finite(padj),
          padj < 0.05
        )
      
      if (nrow(tmp) == 0) next
      
      pos <- tmp %>%
        dplyr::filter(NES > 0) %>%
        dplyr::group_by(annotation, celltype, database) %>%
        dplyr::arrange(dplyr::desc(NES), padj, .by_group = TRUE) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Positive")
      
      neg <- tmp %>%
        dplyr::filter(NES < 0) %>%
        dplyr::group_by(annotation, celltype, database) %>%
        dplyr::arrange(NES, padj, .by_group = TRUE) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Negative")
      
      picked <- dplyr::bind_rows(pos, neg) %>%
        dplyr::ungroup() %>%
        dplyr::select(annotation, celltype, database, direction, pathway, NES, padj, source_file, source_path)
      
      if (nrow(picked) > 0) {
        out_list[[paste0(ann_name, "::", f)]] <- picked
      }
    }
  }
  
  if (length(out_list) == 0) return(data.frame())
  
  dplyr::bind_rows(out_list) %>%
    dplyr::group_by(annotation, celltype, database, direction, pathway) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(annotation, celltype, database, direction) %>%
    dplyr::arrange(
      dplyr::if_else(direction == "Positive", -NES, NES),
      padj,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
}

export_spotcheck_workbook <- function(base_dirs, contrast_filter, out_dir, top_n = 10) {
  inventory <- inspect_result_files(base_dirs, out_dir)
  
  deg_spot  <- collect_detailed_deg_tables(base_dirs, contrast_filter, top_n = top_n)
  gsea_spot <- collect_detailed_gsea_tables(base_dirs, contrast_filter, top_n = top_n)
  
  out_xlsx <- file.path(out_dir, paste0("spotcheck_top", top_n, "_per_celltype.xlsx"))
  
  wb <- openxlsx::createWorkbook()
  
  if (nrow(deg_spot) > 0) {
    openxlsx::addWorksheet(wb, "Top_DEGs")
    openxlsx::writeData(wb, "Top_DEGs", deg_spot)
    openxlsx::freezePane(wb, "Top_DEGs", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Top_DEGs", cols = 1:ncol(deg_spot), widths = "auto")
  }
  
  if (nrow(gsea_spot) > 0) {
    openxlsx::addWorksheet(wb, "Top_Pathways")
    openxlsx::writeData(wb, "Top_Pathways", gsea_spot)
    openxlsx::freezePane(wb, "Top_Pathways", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Top_Pathways", cols = 1:ncol(gsea_spot), widths = "auto")
  }
  
  openxlsx::addWorksheet(wb, "Inventory")
  if (nrow(inventory) > 0) {
    openxlsx::writeData(wb, "Inventory", inventory)
    openxlsx::freezePane(wb, "Inventory", firstRow = TRUE)
  } else {
    openxlsx::writeData(
      wb,
      "Inventory",
      data.frame(message = "No CSV files were found under the base directories.")
    )
  }
  
  if (nrow(deg_spot) == 0 && nrow(gsea_spot) == 0) {
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(
      wb,
      "README",
      data.frame(
        message = c(
          "No significant detailed DEG or pathway results were found by the exporter.",
          "Check the Inventory sheet and result_file_inventory.csv to see what files and columns were detected."
        )
      )
    )
  }
  
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message("Saved spot-check workbook: ", out_xlsx)
  
  invisible(list(
    inventory = inventory,
    deg = deg_spot,
    gsea = gsea_spot,
    out_file = out_xlsx
  ))
}

# =============================================================================
# RUN EXPORT
# =============================================================================
message("Exporting top 10 positive/negative DEGs and pathways per celltype/cluster...")
spotcheck <- export_spotcheck_workbook(
  base_dirs = BASE_DIRS,
  contrast_filter = PAIRED_CONTRAST,
  out_dir = OUT_DIR,
  top_n = 10
)

message("DEG rows exported: ", nrow(spotcheck$deg))
message("GSEA rows exported: ", nrow(spotcheck$gsea))
message("Workbook: ", spotcheck$out_file)



#############GSEA


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(openxlsx)
})

# =============================================================================
# SETTINGS
# =============================================================================
BASE_DIRS <- list(
  Monaco_main  = "deg_gsea_pooled_vs_within_pool_Monaco_main_pruned",
  Monaco_fine  = "deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned",
  Cluster      = "deg_gsea_pooled_vs_within_pool_cluster"
)

OUT_DIR <- "paired_comparison_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAIRED_CONTRAST <- "MildTreated_vs_MildUntreated_paired"

# =============================================================================
# HELPERS
# =============================================================================
safe_read_csv <- function(path) {
  tryCatch(
    data.table::fread(path),
    error = function(e) NULL
  )
}

first_existing_col <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

extract_database_from_filename <- function(file_name) {
  if (grepl("GO_BP", file_name, fixed = TRUE)) return("GO_BP")
  if (grepl("Hallmark", file_name, fixed = TRUE)) return("Hallmark")
  if (grepl("Reactome", file_name, fixed = TRUE)) return("Reactome")
  return("Unknown")
}

extract_celltype_from_subdir <- function(path, anchor_dir) {
  parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/")[[1]]
  idx <- which(parts == anchor_dir)
  if (length(idx) > 0 && idx[length(idx)] < length(parts)) {
    return(parts[idx[length(idx)] + 1])
  }
  return(NA_character_)
}

collect_paired_deg_tables <- function(base_dirs, contrast_filter, top_n = 10) {
  out_list <- list()
  
  for (ann_name in names(base_dirs)) {
    deg_root <- file.path(base_dirs[[ann_name]], "paired_cross_pool", "deg_by_celltype")
    if (!dir.exists(deg_root)) next
    
    csvs <- list.files(
      deg_root,
      pattern = "\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(csvs) == 0) next
    
    for (f in csvs) {
      file_name <- basename(f)
      
      # only keep files for the requested contrast
      if (!grepl(contrast_filter, file_name, fixed = TRUE)) next
      
      df <- safe_read_csv(f)
      if (is.null(df) || nrow(df) == 0) next
      
      gene_col  <- first_existing_col(df, c("gene", "Gene", "symbol", "SYMBOL", "feature", "rowname"))
      logfc_col <- first_existing_col(df, c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC"))
      padj_col  <- first_existing_col(df, c("padj", "FDR", "adj.P.Val", "p_val_adj"))
      
      if (is.na(gene_col) || is.na(logfc_col) || is.na(padj_col)) next
      
      celltype_from_path <- extract_celltype_from_subdir(f, "deg_by_celltype")
      
      tmp <- df %>%
        dplyr::mutate(
          annotation  = ann_name,
          celltype    = celltype_from_path,
          gene        = as.character(.data[[gene_col]]),
          logFC       = suppressWarnings(as.numeric(.data[[logfc_col]])),
          padj        = suppressWarnings(as.numeric(.data[[padj_col]])),
          source_file = file_name,
          source_path = f
        ) %>%
        dplyr::filter(
          !is.na(celltype), celltype != "",
          !is.na(gene), gene != "",
          !is.na(logFC), is.finite(logFC),
          !is.na(padj), is.finite(padj),
          padj < 0.05
        )
      
      if (nrow(tmp) == 0) next
      
      pos <- tmp %>%
        dplyr::filter(logFC > 0) %>%
        dplyr::arrange(dplyr::desc(logFC), padj) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Positive")
      
      neg <- tmp %>%
        dplyr::filter(logFC < 0) %>%
        dplyr::arrange(logFC, padj) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Negative")
      
      picked <- dplyr::bind_rows(pos, neg) %>%
        dplyr::select(
          annotation, celltype, direction,
          gene, logFC, padj,
          source_file, source_path
        )
      
      if (nrow(picked) > 0) {
        out_list[[paste0(ann_name, "::", f)]] <- picked
      }
    }
  }
  
  if (length(out_list) == 0) return(data.frame())
  
  dplyr::bind_rows(out_list) %>%
    dplyr::group_by(annotation, celltype, direction, gene) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(annotation, celltype, direction) %>%
    dplyr::arrange(
      dplyr::if_else(direction == "Positive", -logFC, logFC),
      padj,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
}

collect_paired_gsea_tables <- function(base_dirs, contrast_filter, top_n = 10) {
  out_list <- list()
  
  for (ann_name in names(base_dirs)) {
    gsea_root <- file.path(base_dirs[[ann_name]], "paired_cross_pool", "gsea_by_celltype")
    if (!dir.exists(gsea_root)) next
    
    csvs <- list.files(
      gsea_root,
      pattern = "_results(_sig)?\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(csvs) == 0) next
    
    for (f in csvs) {
      file_name <- basename(f)
      
      # only keep files for the requested contrast
      if (!grepl(contrast_filter, file_name, fixed = TRUE)) next
      
      df <- safe_read_csv(f)
      if (is.null(df) || nrow(df) == 0) next
      
      pathway_col <- first_existing_col(df, c("Description", "pathway", "Pathway", "term", "ID"))
      nes_col     <- first_existing_col(df, c("NES", "nes"))
      padj_col    <- first_existing_col(df, c("padj", "p.adjust", "FDR", "qvalue"))
      
      if (is.na(pathway_col) || is.na(nes_col) || is.na(padj_col)) next
      
      celltype_from_path <- extract_celltype_from_subdir(f, "gsea_by_celltype")
      db_from_file <- extract_database_from_filename(file_name)
      
      tmp <- df %>%
        dplyr::mutate(
          annotation  = ann_name,
          celltype    = celltype_from_path,
          database    = db_from_file,
          pathway     = as.character(.data[[pathway_col]]),
          NES         = suppressWarnings(as.numeric(.data[[nes_col]])),
          padj        = suppressWarnings(as.numeric(.data[[padj_col]])),
          source_file = file_name,
          source_path = f
        ) %>%
        dplyr::filter(
          !is.na(celltype), celltype != "",
          !is.na(database), database != "",
          !is.na(pathway), pathway != "",
          !is.na(NES), is.finite(NES),
          !is.na(padj), is.finite(padj),
          padj < 0.05
        )
      
      if (nrow(tmp) == 0) next
      
      pos <- tmp %>%
        dplyr::filter(NES > 0) %>%
        dplyr::arrange(dplyr::desc(NES), padj) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Positive")
      
      neg <- tmp %>%
        dplyr::filter(NES < 0) %>%
        dplyr::arrange(NES, padj) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(direction = "Negative")
      
      picked <- dplyr::bind_rows(pos, neg) %>%
        dplyr::select(
          annotation, celltype, database, direction,
          pathway, NES, padj,
          source_file, source_path
        )
      
      if (nrow(picked) > 0) {
        out_list[[paste0(ann_name, "::", f)]] <- picked
      }
    }
  }
  
  if (length(out_list) == 0) return(data.frame())
  
  dplyr::bind_rows(out_list) %>%
    dplyr::group_by(annotation, celltype, database, direction, pathway) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(annotation, celltype, database, direction) %>%
    dplyr::arrange(
      dplyr::if_else(direction == "Positive", -NES, NES),
      padj,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
}

export_spotcheck_workbook <- function(base_dirs, contrast_filter, out_dir, top_n = 10) {
  deg_spot  <- collect_paired_deg_tables(base_dirs, contrast_filter, top_n = top_n)
  gsea_spot <- collect_paired_gsea_tables(base_dirs, contrast_filter, top_n = top_n)
  
  out_xlsx <- file.path(
    out_dir,
    paste0("spotcheck_top", top_n, "_", contrast_filter, ".xlsx")
  )
  
  wb <- openxlsx::createWorkbook()
  
  if (nrow(deg_spot) > 0) {
    openxlsx::addWorksheet(wb, "Top_DEGs")
    openxlsx::writeData(wb, "Top_DEGs", deg_spot)
    openxlsx::freezePane(wb, "Top_DEGs", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Top_DEGs", cols = 1:ncol(deg_spot), widths = "auto")
  }
  
  if (nrow(gsea_spot) > 0) {
    openxlsx::addWorksheet(wb, "Top_Pathways")
    openxlsx::writeData(wb, "Top_Pathways", gsea_spot)
    openxlsx::freezePane(wb, "Top_Pathways", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Top_Pathways", cols = 1:ncol(gsea_spot), widths = "auto")
  }
  
  if (nrow(deg_spot) == 0 && nrow(gsea_spot) == 0) {
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(
      wb,
      "README",
      data.frame(
        message = paste("No paired DEG or GSEA results found for contrast:", contrast_filter)
      )
    )
  }
  
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message("Saved: ", out_xlsx)
  
  invisible(list(
    deg = deg_spot,
    gsea = gsea_spot,
    out_file = out_xlsx
  ))
}

# =============================================================================
# RUN
# =============================================================================
spotcheck <- export_spotcheck_workbook(
  base_dirs = BASE_DIRS,
  contrast_filter = PAIRED_CONTRAST,
  out_dir = OUT_DIR,
  top_n = 10
)

message("DEG rows exported: ", nrow(spotcheck$deg))
message("GSEA rows exported: ", nrow(spotcheck$gsea))
message("Workbook: ", spotcheck$out_file)