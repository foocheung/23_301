# ============================================================================
# REUSABLE QMD HELPER FUNCTIONS FOR SINGLE-CELL ANALYSIS
# ============================================================================
# These functions create tabsets that automatically iterate through metadata
# columns to generate UMAPs, violin plots, heatmaps, and tables.
# ============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(DT)
library(scales)

# ----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ----------------------------------------------------------------------------

#' Get plottable metadata columns
#' 
#' @param obj Seurat object
#' @param max_levels Maximum number of unique values for categorical variables
#' @param exclude_cols Character vector of column names to exclude
#' @return Character vector of plottable column names
get_plottable_columns <- function(obj, max_levels = 50, exclude_cols = NULL) {
  meta <- obj@meta.data
  cols <- colnames(meta)
  
  if (!is.null(exclude_cols)) {
    cols <- setdiff(cols, exclude_cols)
  }
  
  plottable <- vapply(cols, function(cn) {
    v <- meta[[cn]]
    is.atomic(v) && !is.list(v) &&
      length(na.omit(unique(v))) > 1L &&
      length(na.omit(unique(v))) <= max_levels
  }, logical(1))
  
  names(plottable)[plottable]
}

#' Convert to discrete factor for plotting
#' 
#' @param x Vector to convert
#' @param n_bins Number of bins for numeric data
#' @return Factor
as_discrete_factor <- function(x, n_bins = 5) {
  if (is.factor(x) || is.character(x) || is.logical(x)) return(factor(x))
  qs <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(qs) < 3) return(factor(ifelse(is.na(x), NA, "value")))
  cut(x, breaks = qs, include.lowest = TRUE, dig.lab = 6)
}

#' Get discrete color palette
#' 
#' @param n Number of colors needed
#' @return Character vector of hex colors
discrete_palette <- function(n) {
  if (n == 1) return("#4DBBD5")
  if (n == 2) return(c("#E64B35", "#4DBBD5"))
  scales::hue_pal()(n)
}

# ----------------------------------------------------------------------------
# UMAP PLOTTING HELPERS (fixed API)
# ----------------------------------------------------------------------------

# palette + factor helpers (keep if you already have them)
discrete_palette <- function(n) {
  if (n <= 1) return("#4DBBD5")
  if (n == 2) return(c("#E64B35", "#4DBBD5"))
  scales::hue_pal()(n)
}

as_discrete_factor <- function(x, n_bins = 5) {
  if (is.factor(x) || is.character(x) || is.logical(x)) return(factor(x))
  qs <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(qs) < 3) return(factor(ifelse(is.na(x), NA, "value")))
  cut(x, breaks = qs, include.lowest = TRUE, dig.lab = 6)
}

# Small wrapper to use scCustomize when present, else Seurat::DimPlot
.dimplot_umap <- function(obj, group.by, split.by = NULL, pt.size = 0.3,
                          colors_use = NULL, label = FALSE, num_columns = 2,
                          reduction = "umap") {
  use_scCustom <- exists("DimPlot_scCustom")
  if (use_scCustom) {
    if (is.null(split.by)) {
      p <- DimPlot_scCustom(
        obj, reduction = reduction,
        group.by = group.by,
        pt.size = pt.size,
        label = label,
        colors_use = colors_use
      )
    } else {
      p <- DimPlot_scCustom(
        obj, reduction = reduction,
        split.by = split.by,
        group.by = group.by,
        pt.size = pt.size,
        label = label,
        num_columns = num_columns,
        colors_use = colors_use
      )
    }
  } else {
    # Seurat::DimPlot uses `cols` instead of `colors_use`
    cols <- colors_use
    if (is.null(split.by)) {
      p <- Seurat::DimPlot(
        obj, reduction = reduction,
        group.by = group.by,
        pt.size = pt.size,
        label = label,
        cols = cols
      )
    } else {
      p <- Seurat::DimPlot(
        obj, reduction = reduction,
        split.by = split.by,
        group.by = group.by,
        pt.size = pt.size,
        label = label,
        ncol = num_columns,
        cols = cols
      )
    }
  }
  p
}

# ----------------------------------------------------------------------------
# UMAP PLOTTING FUNCTIONS (no duplicates; unified signatures)
# ----------------------------------------------------------------------------

# UMAP colored by a single metadata column (no split)
plot_umap_by_meta <- function(obj, colname, pt.size = 0.8, reduction = "umap",
                              colors_use = NULL, label = FALSE) {
  v <- obj@meta.data[[colname]]
  if (all(is.na(v)) || length(na.omit(unique(v))) <= 1L) return(NULL)
  
  fac <- as_discrete_factor(v)
  nlev <- length(levels(fac)); if (nlev <= 1L) return(NULL)
  
  cols <- if (!is.null(colors_use)) colors_use else discrete_palette(nlev)
  
  tmp_col <- paste0(".TMP_", colname)
  obj@meta.data[[tmp_col]] <- fac
  
  .dimplot_umap(
    obj, group.by = tmp_col, split.by = NULL,
    pt.size = pt.size, colors_use = cols, label = label, reduction = reduction
  ) +
    labs(title = paste0("UMAP by ", colname), color = colname) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

# UMAP with fixed split.by and variable group.by
plot_umap_by_meta_split <- function(obj, group_col, split_col,
                                    pt.size = 0.8, reduction = "umap",
                                    colors_use = NULL, num_columns = 2,
                                    label = TRUE) {
  v <- obj@meta.data[[group_col]]
  if (all(is.na(v)) || length(na.omit(unique(v))) <= 1L) return(NULL)
  
  fac <- as_discrete_factor(v)
  nlev <- length(levels(fac)); if (nlev <= 1L) return(NULL)
  
  cols <- if (!is.null(colors_use)) colors_use else discrete_palette(nlev)
  
  tmp_col <- paste0(".TMP_", group_col)
  obj@meta.data[[tmp_col]] <- fac
  
  .dimplot_umap(
    obj, group.by = tmp_col, split.by = split_col,
    pt.size = pt.size, colors_use = cols, label = label,
    num_columns = num_columns, reduction = reduction
  ) +
    labs(title = paste0("UMAP (split by ", split_col, ") — colored by ", group_col),
         color = group_col) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}

# Tabset of UMAPs for all plottable metadata columns
generate_umap_tabset <- function(obj, exclude_cols = NULL, 
                                 priority_cols = NULL, split_by = NULL,
                                 pt.size = 0.8, colors_use = NULL,
                                 reduction = "umap", num_columns = 2,
                                 label = TRUE, max_tabs = NULL) {
  plottable <- get_plottable_columns(obj, exclude_cols = exclude_cols)
  
  if (!is.null(priority_cols)) {
    priority_cols <- intersect(priority_cols, plottable)
    other_cols <- setdiff(plottable, priority_cols)
    plottable <- c(priority_cols, other_cols)
  }
  
  # limit for testing
  if (!is.null(max_tabs)) {
    plottable <- head(plottable, max_tabs)
  }
  
  cat("\n::: {.panel-tabset}\n\n")
  for (col in plottable) {
    cat("### ", col, "\n\n", sep = "")
    
    p <- if (is.null(split_by)) {
      plot_umap_by_meta(
        obj, col, pt.size = pt.size, reduction = reduction,
        colors_use = colors_use, label = label
      )
    } else {
      plot_umap_by_meta_split(
        obj, col, split_col = split_by, pt.size = pt.size,
        reduction = reduction, colors_use = colors_use,
        num_columns = num_columns, label = label
      )
    }
    
    if (is.null(p)) {
      cat("*No variation or not plottable*\n\n")
    } else {
      print(p)
      cat("\n\n")
    }
  }
  cat(":::\n\n")
}


# ----------------------------------------------------------------------------
# VIOLIN PLOT FUNCTIONS
# ----------------------------------------------------------------------------

#' Generate violin plot for a feature colored by metadata
#' 
#' @param obj Seurat object
#' @param features Features to plot (genes or QC metrics)
#' @param group_by Metadata column for grouping
#' @param split_by Optional metadata column for splitting
#' @param pt.size Point size
#' @param log_scale Use log scale for y-axis
#' @return ggplot object
plot_violin_by_meta <- function(obj, features, group_by, split_by = NULL,
                                pt.size = 0, log_scale = FALSE) {
  p <- VlnPlot(obj, features = features, group.by = group_by,
               split.by = split_by, pt.size = pt.size)
  
  if (log_scale) {
    p <- p + scale_y_log10(labels = comma)
  }
  
  p + 
    labs(title = paste0(features, " by ", group_by),
         x = NULL, y = features) +
    theme_minimal(base_size = 12) +
    theme(legend.position = if(is.null(split_by)) "none" else "bottom")
}


# ---------- Helpers for palettes ----------
discrete_palette <- function(n) {
  if (n == 1) return("#4DBBD5")
  if (n == 2) return(c("#E64B35", "#4DBBD5"))
  scales::hue_pal()(n)
}

# ---------- Core violin plot ----------
# Plots a single feature as a violin, grouped by a metadata column, optionally split.
plot_violin_by_meta <- function(obj, feature, group_col, split_by = NULL, log_scale = FALSE) {
  # sanity checks
  if (!feature %in% rownames(obj)) stop(sprintf("Feature '%s' not found in object.", feature))
  if (!group_col %in% colnames(obj@meta.data)) stop(sprintf("group_col '%s' not in metadata.", group_col))
  if (!is.null(split_by) && !split_by %in% colnames(obj@meta.data)) {
    stop(sprintf("split_by '%s' not in metadata.", split_by))
  }
  
  # determine palette from number of levels in group_col
  levs <- unique(na.omit(obj@meta.data[[group_col]]))
  levs <- levs[order(as.character(levs))]
  nlev <- length(levs)
  if (nlev <= 1L) stop(sprintf("Grouping column '%s' has only one level.", group_col))
  cols <- discrete_palette(nlev)
  
  p <- VlnPlot(
    object   = obj,
    features = feature,
    group.by = group_col,
    split.by = split_by,
    pt.size  = 0,
    cols     = cols
  ) +
    labs(
      title = sprintf("Violin: %s (grouped by %s%s)",
                      feature, group_col,
                      if (!is.null(split_by)) paste0(", split by ", split_by) else ""),
      x = NULL, y = "Expression"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  if (isTRUE(log_scale)) {
    # log10 scale; keep zeros visible by adding +1 pseudocount in the label only
    p <- p + scale_y_log10() +
      labs(subtitle = "Y-axis on log10 scale")
  }
  
  p
}

# ---------- Nested tabset generator ----------
#' Generate nested tabsets of violins.
#' Outer tabs = group_by_cols; inner tabs = features (each with its own plot).
#' Keeps split_by fixed across all plots if provided.
#'
#' @param obj Seurat object
#' @param features character vector of gene features (rows of obj)
#' @param group_by_cols character vector of metadata columns to group by (outer tabs)
#' @param split_by single metadata column to split violins by (optional, fixed)
#' @param log_scale named logical vector (e.g., c(GENE1=TRUE, GENE2=FALSE)); TRUE => log10 y
#' @return prints markdown tabsets
# ---------- palettes ----------
discrete_palette <- function(n) {
  if (n == 1) return("#4DBBD5")
  if (n == 2) return(c("#E64B35", "#4DBBD5"))
  scales::hue_pal()(n)
}

# ---------- core violin plot (genes OR metadata) ----------
plot_violin_by_meta <- function(obj, feature, group_col, split_by = NULL, log_scale = FALSE) {
  # sanity checks for grouping/splitting columns
  stopifnot(group_col %in% colnames(obj@meta.data))
  if (!is.null(split_by)) stopifnot(split_by %in% colnames(obj@meta.data))
  
  # Determine number of groups for color palette (x-axis groups)
  levs <- unique(na.omit(obj@meta.data[[group_col]]))
  levs <- levs[order(as.character(levs))]
  nlev <- length(levs)
  if (nlev <= 1L) stop(sprintf("Grouping column '%s' has only one level.", group_col))
  cols <- discrete_palette(nlev)
  
  # Try Seurat's VlnPlot directly (it supports metadata names as features, e.g., nCount_RNA)
  p <- tryCatch({
    Seurat::VlnPlot(
      object   = obj,
      features = feature,    # works for both gene features and metadata columns
      group.by = group_col,
      split.by = split_by,
      pt.size  = 0,
      cols     = cols
    )
  }, error = function(e) NULL)
  
  # Fallback: manual ggplot if VlnPlot fails for some reason
  if (is.null(p)) {
    # fetch y-values: gene expression or metadata
    yvals <- if (feature %in% rownames(obj)) {
      as.numeric(Seurat::GetAssayData(obj, slot = "data")[feature, , drop = TRUE])
    } else if (feature %in% colnames(obj@meta.data)) {
      obj@meta.data[[feature]]
    } else {
      stop(sprintf("Feature '%s' not found as gene or metadata.", feature))
    }
    if (!is.numeric(yvals)) {
      stop(sprintf("Feature '%s' is not numeric; cannot plot violin.", feature))
    }
    
    df <- data.frame(
      value   = yvals,
      group   = obj@meta.data[[group_col]],
      splitby = if (!is.null(split_by)) obj@meta.data[[split_by]] else factor("all")
    )
    df <- df[!is.na(df$value) & !is.na(df$group) & !is.na(df$splitby), , drop = FALSE]
    
    p <- ggplot(df, aes(x = group, y = value, fill = if (!is.null(split_by)) splitby else group)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.8, position = position_dodge(width = 0.9)) +
      scale_fill_manual(values = if (is.null(split_by)) cols else discrete_palette(length(unique(df$splitby)))) +
      labs(
        title = sprintf("Violin: %s (grouped by %s%s)",
                        feature, group_col,
                        if (!is.null(split_by)) paste0(", split by ", split_by) else ""),
        x = NULL, y = "Expression"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
  }
  
  if (isTRUE(log_scale)) {
    # Note: log10 requires positive values
    p <- p + scale_y_log10() +
      labs(subtitle = "Y-axis on log10 scale")
  }
  
  p
}

# ---------- nested tabset generator ----------
generate_violin_tabset <- function(obj, features, group_by_cols,
                                   split_by = NULL, log_scale = NULL) {
  cat("\n:::{.panel-tabset}\n\n")  # outer tabset: one tab per group_col
  
  for (group_col in group_by_cols) {
    cat("### ", group_col, "\n\n", sep = "")
    
    # inner tabset: one tab per feature
    cat(":::{.panel-tabset}\n\n")
    
    for (feat in features) {
      # safe lookup: returns NA if missing; isTRUE(NA) => FALSE
      use_log <- isTRUE(unname(log_scale[feat])[1])
      
      cat("#### ", feat, "\n\n", sep = "")
      
      p <- tryCatch(
        plot_violin_by_meta(obj, feat, group_col, split_by = split_by, log_scale = use_log),
        error = function(e) {
          cat(sprintf("**Skipped**: %s\n\n", e$message))
          NULL
        }
      )
      if (!is.null(p)) {
        print(p)
        cat("\n\n")
      }
    }
    
    cat(":::\n\n")  # close inner tabset
  }
  
  cat(":::\n\n")    # close outer tabset
}


#' #' Generate tabset of violin plots for multiple features
#' #' 
#' #' @param obj Seurat object
#' #' @param features Vector of features to plot
#' #' @param group_by_cols Metadata columns to group by
#' #' @param split_by Optional split variable
#' #' @param log_scale Named vector of features requiring log scale
#' #' @return Prints markdown tabset
#' generate_violin_tabset <- function(obj, features, group_by_cols, 
#'                                    split_by = NULL, log_scale = NULL) {
#'   cat("\n::: {.panel-tabset}\n\n")
#'   
#'   for (group_col in group_by_cols) {
#'     cat("###", group_col, "\n\n")
#'     
#'     for (feat in features) {
#'       use_log <- !is.null(log_scale) && feat %in% names(log_scale) && log_scale[feat]
#'       
#'       p <- tryCatch({
#'         plot_violin_by_meta(obj, feat, group_col, split_by, log_scale = use_log)
#'       }, error = function(e) NULL)
#'       
#'       if (!is.null(p)) {
#'         print(p)
#'         cat("\n\n")
#'       }
#'     }
#'   }
#'   
#'   cat(":::\n\n")
#' }

# ----------------------------------------------------------------------------
# HEATMAP FUNCTIONS
# ----------------------------------------------------------------------------

#' Generate count heatmap for two categorical variables
#' 
#' @param obj Seurat object
#' @param x_var X-axis variable
#' @param y_var Y-axis variable
#' @param facet_var Optional faceting variable
#' @param log_scale Use log scale for fill
#' @return ggplot object
plot_count_heatmap <- function(obj, x_var, y_var, facet_var = NULL, 
                               log_scale = TRUE) {
  count_data <- if (is.null(facet_var)) {
    table(obj@meta.data[[y_var]], obj@meta.data[[x_var]])
  } else {
    table(obj@meta.data[[y_var]], 
          obj@meta.data[[x_var]], 
          obj@meta.data[[facet_var]])
  }
  
  df <- as.data.frame(count_data) %>%
    dplyr::filter(Freq > 0)
  
  p <- ggplot(df, aes(x = .data[[names(df)[2]]], 
                      y = .data[[names(df)[1]]], 
                      fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Freq), size = 3, color = "black")
  
  if (log_scale) {
    p <- p + scale_fill_gradient(low = "white", high = "steelblue", 
                                 trans = "log10", labels = comma)
  } else {
    p <- p + scale_fill_gradient(low = "white", high = "steelblue", 
                                 labels = comma)
  }
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", names(df)[3])), scales = "free")
  }
  
  p + 
    labs(x = x_var, y = y_var, fill = "Count") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

#' Generate expression heatmap
#' 
#' @param obj Seurat object
#' @param features Genes to plot
#' @param group_by Metadata column for grouping
#' @param slot Data slot to use
#' @param scale_rows Scale by row (z-score)
#' @return ggplot object
plot_expression_heatmap <- function(obj, features, group_by, 
                                    slot = "data", scale_rows = TRUE) {
  # Get available features
  features <- intersect(features, rownames(obj))
  if (length(features) == 0) return(NULL)
  
  # Average expression per group
  avg_expr <- AverageExpression(obj, features = features,
                                group.by = group_by,
                                slot = slot)[[DefaultAssay(obj)]]
  
  
  n_na <- sum(is.na(avg_expr))
  n_inf <- sum(is.infinite(avg_expr))
  n_nan <- sum(is.nan(avg_expr))
  
  if (n_na > 0 || n_inf > 0 || n_nan > 0) {
    log_warn("Found problematic values: {n_na} NAs, {n_inf} Inf, {n_nan} NaN")
    
    # Replace NA/NaN with 0
    avg_expr[is.na(avg_expr) | is.nan(avg_expr)] <- 0
    
    # Replace Inf with max finite value
    if (any(is.infinite(avg_expr))) {
      max_finite <- max(avg_expr[is.finite(avg_expr)], na.rm = TRUE)
      avg_expr[is.infinite(avg_expr) & avg_expr > 0] <- max_finite
      avg_expr[is.infinite(avg_expr) & avg_expr < 0] <- -max_finite
    }
  }
  
  # Remove genes with zero variance (constant across all groups)
  gene_var <- apply(avg_expr, 1, var, na.rm = TRUE)
  zero_var_genes <- names(gene_var)[gene_var == 0 | is.na(gene_var)]
  
  if (length(zero_var_genes) > 0) {
    log_warn("Removing {length(zero_var_genes)} genes with zero variance")
    avg_expr <- avg_expr[!rownames(avg_expr) %in% zero_var_genes, , drop = FALSE]
  }
  
  # Final check
  if (nrow(avg_expr) == 0) {
    stop("No genes remaining after filtering")
  }
  
  if (any(!is.finite(avg_expr))) {
    log_error("Still have non-finite values after cleaning!")
    avg_expr[!is.finite(avg_expr)] <- 0
  }
  
  # Convert to long format
  df <- as.data.frame(avg_expr) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "group", values_to = "expression")
  
  if (scale_rows) {
    df <- df %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(scaled_expr = scale(expression)[,1]) %>%
      dplyr::ungroup()
    fill_var <- "scaled_expr"
    fill_name <- "Scaled\nExpression"
  } else {
    fill_var <- "expression"
    fill_name <- "Expression"
  }
  
  ggplot(df, aes(x = group, y = gene, fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                         midpoint = 0, name = fill_name) +
    labs(x = group_by, y = "Gene") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

# ----------------------------------------------------------------------------
# TABLE FUNCTIONS
# ----------------------------------------------------------------------------

#' Create summary table from metadata
#' 
#' @param obj Seurat object
#' @param group_vars Variables to group by
#' @param summary_vars Variables to summarize
#' @param summary_funs Named list of summary functions
#' @return tibble
create_summary_table <- function(obj, group_vars, summary_vars = NULL, 
                                 summary_funs = NULL) {
  df <- obj@meta.data
  
  if (is.null(summary_vars)) {
    # Just count
    df %>%
      dplyr::group_by(across(all_of(group_vars))) %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
      dplyr::arrange(desc(n))
  } else {
    # Custom summaries
    if (is.null(summary_funs)) {
      summary_funs <- list(
        median = ~median(.x, na.rm = TRUE),
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      )
    }
    
    df %>%
      dplyr::group_by(across(all_of(group_vars))) %>%
      dplyr::summarise(
        n = n(),
        across(all_of(summary_vars), summary_funs, .names = "{.col}_{.fn}"),
        .groups = "drop"
      )
  }
}

#' Render interactive datatable
#' 
#' @param data Data frame or tibble
#' @param caption Table caption
#' @param page_length Rows per page
#' @param buttons Export buttons to show
#' @param filename Base filename for exports
#' @return DT datatable
render_datatable <- function(data, caption = NULL, page_length = 25,
                             buttons = c("copy", "csv", "excel"),
                             filename = "data_export") {
  datatable(
    data,
    extensions = c("Buttons", "Scroller"),
    options = list(
      dom = "Bfrtip",
      buttons = lapply(buttons, function(btn) {
        if (btn %in% c("csv", "excel")) {
          list(extend = btn, filename = filename)
        } else {
          btn
        }
      }),
      scrollX = TRUE,
      deferRender = TRUE,
      scrollY = 500,
      scroller = TRUE,
      pageLength = page_length
    ),
    caption = caption,
    filter = "top",
    rownames = FALSE
  )
}

# ----------------------------------------------------------------------------
# PROPORTION ANALYSIS
# ----------------------------------------------------------------------------

#' Calculate proportions across groups
#' 
#' @param obj Seurat object
#' @param cluster_col Cell type column
#' @param subject_col Subject identifier
#' @param condition_col Condition/treatment column
#' @return tibble with proportions
calculate_proportions <- function(obj, cluster_col, subject_col, condition_col) {
  obj@meta.data %>%
    dplyr::count(.data[[subject_col]], .data[[condition_col]], .data[[cluster_col]]) %>%
    dplyr::group_by(.data[[subject_col]], .data[[condition_col]]) %>%
    dplyr::mutate(
      total = sum(n),
      proportion = n / total
    ) %>%
    dplyr::ungroup()
}

#' Plot paired proportion changes
#' 
#' @param prop_data Output from calculate_proportions()
#' @param cluster_col Cell type column name
#' @param subject_col Subject column name
#' @param condition_col Condition column name
#' @param ncol Number of facet columns
#' @return ggplot object
plot_paired_proportions <- function(prop_data, cluster_col, subject_col, 
                                    condition_col, ncol = 4) {
  ggplot(prop_data, 
         aes(x = .data[[condition_col]], 
             y = proportion,
             group = .data[[subject_col]], 
             color = .data[[subject_col]])) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 2.5) +
    facet_wrap(as.formula(paste("~", cluster_col)), 
               scales = "free_y", ncol = ncol) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = "Cell Type Proportion Changes (Paired by Subject)",
         subtitle = "Lines connect the same subject across conditions",
         x = NULL, y = "Proportion within Subject",
         color = subject_col) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey90"),
          panel.grid.minor = element_blank())
}

# ----------------------------------------------------------------------------
# USAGE EXAMPLES (commented out)
# ----------------------------------------------------------------------------

# # In your QMD file, use like this:
# 
# ```{r umap_all_metadata, results='asis', fig.width=10, fig.height=8}
# generate_umap_tabset(
#   obj, 
#   priority_cols = c("Treatment", "SubjectTimePoint_ordered"),
#   exclude_cols = c("orig.ident", "nCount_RNA", "nFeature_RNA")
# )
# ```
# 
# ```{r umap_split_metadata, results='asis', fig.width=10, fig.height=12}
# generate_umap_tabset(
#   obj,
#   split_by = "SubjectTimePoint_ordered",
#   priority_cols = c("Treatment", "PanGI_L2")
# )
# ```
# 
# ```{r violin_qc, results='asis', fig.width=12, fig.height=6}
# generate_violin_tabset(
#   obj,
#   features = c("nCount_RNA", "nFeature_RNA", "percent.mito"),
#   group_by_cols = c("Treatment", "Subject"),
#   split_by = "Treatment",
#   log_scale = c(nCount_RNA = TRUE, nFeature_RNA = TRUE)
# )
# ```
# 
# ```{r summary_table}
# summary_tbl <- create_summary_table(
#   obj,
#   group_vars = c("Subject", "Treatment"),
#   summary_vars = c("nCount_RNA", "nFeature_RNA", "percent.mito")
# )
# 
# render_datatable(
#   summary_tbl,
#   caption = "Sample-level QC summary",
#   filename = "qc_summary"
# )
# ```