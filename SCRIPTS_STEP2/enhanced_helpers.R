# ============================================================================
# Enhanced Helper Functions for Single-Cell RNA-seq Analysis
# ============================================================================

# Load required libraries
#suppressPackageStartageMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(Seurat)
  library(logger)
  library(progress)
#})

# ============================================================================
# 1. CONFIGURATION & VALIDATION
# ============================================================================

#' Validate configuration file
#' @param config List from yaml::read_yaml()
#' @return NULL (throws errors if invalid)
`%R%` <- function(char, times) {
  strrep(char, times)
}

validate_config <- function(config) {
  log_info("Validating configuration...")
  
  # Check input file exists
  if (!file.exists(config$input$rds_file)) {
    log_error("Input RDS file not found: {config$input$rds_file}")
    stop("Input file does not exist")
  }
  
  # Validate QC thresholds
  if (config$qc$nfeature_min >= config$qc$nfeature_max) {
    stop("nfeature_min must be less than nfeature_max")
  }
  if (config$qc$ncount_min >= config$qc$ncount_max) {
    stop("ncount_min must be less than ncount_max")
  }
  
  # Validate FDR cutoff
  if (config$analysis$fdr_cutoff <= 0 || config$analysis$fdr_cutoff >= 1) {
    stop("FDR cutoff must be between 0 and 1")
  }
  
  # Check parallel cores
#  if (config$advanced$parallel_cores > parallel::detectCores()) {
#    log_warn("Requested cores ({config$advanced$parallel_cores}) exceeds available ({parallel::detectCores()})")
#  }
  
  log_success("Configuration validated successfully")
  invisible(NULL)
}

#' Setup logging with file output
#' @param output_dir Output directory path
#' @param log_level Logging level (INFO, DEBUG, WARN, ERROR)
setup_logging <- function(output_dir, log_level = "INFO") {
  dir.create(file.path(output_dir, "logs"), showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(output_dir, "logs", 
                        sprintf("analysis_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
  
  log_appender(appender_tee(log_file))
  log_threshold(log_level)
  
  log_info("=" %R% 70)
  log_info("Single-Cell RNA-seq Analysis Pipeline")
  log_info("Started: {Sys.time()}")
  log_info("=" %R% 70)
  
  return(log_file)
}

# ============================================================================
# 2. ADVANCED QC FUNCTIONS
# ============================================================================

#' Generate interactive QC scatter plots
#' @param obj Seurat object
#' @param color_by Metadata column to color by
#' @return plotly object

plot_qc_interactive <- function(obj, color_by = "Treatment", annotation_col = NULL) {
  meta <- obj@meta.data
  
  # Build hover text dynamically
  hover_text <- paste0(
    "Cell: ", rownames(meta), "<br>",
    "UMI: ", scales::comma(meta$nCount_RNA), "<br>",
    "Genes: ", scales::comma(meta$nFeature_RNA), "<br>",
    "Mito %: ", round(meta$percent.mito, 2), "<br>"
  )
  
  # Add ribo if available
  if ("percent.ribo" %in% colnames(meta) && !all(is.na(meta$percent.ribo))) {
    hover_text <- paste0(hover_text, "Ribo %: ", round(meta$percent.ribo, 2), "<br>")
  }
  
  # Add novelty if available
  if ("novelty_score" %in% colnames(meta) && !all(is.na(meta$novelty_score))) {
    hover_text <- paste0(hover_text, "Novelty: ", round(meta$novelty_score, 3), "<br>")
  }
  
  # Add annotation if specified and column exists
  if (!is.null(annotation_col) && annotation_col %in% colnames(meta)) {
    hover_text <- paste0(
      hover_text,
      "Cell Type: ", meta[[annotation_col]], "<br>"
    )
  }
  
  # Add color_by variable
  if (color_by %in% colnames(meta)) {
    hover_text <- paste0(
      hover_text,
      color_by, ": ", meta[[color_by]]
    )
  }
  
  p <- plot_ly(
    data = meta,
    x = ~nCount_RNA,
    y = ~nFeature_RNA,
    color = as.formula(paste0("~", color_by)),
    type = "scatter",
    mode = "markers",
    marker = list(size = 3, opacity = 0.6),
    text = hover_text,
    hoverinfo = "text"
  ) %>%
    layout(
      title = "Interactive QC: UMI vs Genes Detected",
      xaxis = list(title = "UMI Counts", type = "log"),
      yaxis = list(title = "Genes Detected", type = "log"),
      hovermode = "closest"
    )
  
  return(p)
}

#' Generate interactive QC plot with multiple annotation columns
#' @param obj Seurat object
#' @param color_by Metadata column to color by
#' @param annotation_cols Vector of annotation columns to include in hover
#' @return plotly object
plot_qc_interactive_multi <- function(obj, 
                                      color_by = "Treatment",
                                      annotation_cols = NULL) {
  meta <- obj@meta.data
  
  # Build hover text
  hover_text <- paste0(
    "Cell: ", rownames(meta), "<br>",
    "UMI: ", scales::comma(meta$nCount_RNA), "<br>",
    "Genes: ", scales::comma(meta$nFeature_RNA), "<br>",
    "Mito %: ", round(meta$percent.mito, 2), "<br>"
  )
  
  # Add ribo if available
  if ("percent.ribo" %in% colnames(meta) && !all(is.na(meta$percent.ribo))) {
    hover_text <- paste0(hover_text, "Ribo %: ", round(meta$percent.ribo, 2), "<br>")
  }
  
  # Add novelty if available
  if ("novelty_score" %in% colnames(meta) && !all(is.na(meta$novelty_score))) {
    hover_text <- paste0(hover_text, "Novelty: ", round(meta$novelty_score, 3), "<br>")
  }
  
  # Add all requested annotation columns
  if (!is.null(annotation_cols)) {
    for (col in annotation_cols) {
      if (col %in% colnames(meta)) {
        col_label <- gsub("_", " ", col)  # Make prettier
        col_label <- tools::toTitleCase(col_label)
        hover_text <- paste0(
          hover_text,
          col_label, ": ", meta[[col]], "<br>"
        )
      }
    }
  }
  
  # Add color_by variable last
  if (color_by %in% colnames(meta)) {
    hover_text <- paste0(
      hover_text,
      color_by, ": ", meta[[color_by]]
    )
  }
  
  p <- plot_ly(
    data = meta,
    x = ~nCount_RNA,
    y = ~nFeature_RNA,
    color = as.formula(paste0("~", color_by)),
    type = "scatter",
    mode = "markers",
    marker = list(size = 3, opacity = 0.6),
    text = hover_text,
    hoverinfo = "text"
  ) %>%
    layout(
      title = "Interactive QC: UMI vs Genes Detected",
      xaxis = list(title = "UMI Counts", type = "log"),
      yaxis = list(title = "Genes Detected", type = "log"),
      hovermode = "closest"
    )
  
  return(p)
}

calculate_complexity <- function(obj) {
  meta <- obj@meta.data
  
  complexity <- meta %>%
    mutate(
      novelty_score = nFeature_RNA / nCount_RNA,
      log10_genes = log10(nFeature_RNA),
      log10_umi = log10(nCount_RNA)
    )
  
  obj@meta.data <- complexity
  return(obj)
}

#' QC failure analysis
#' @param obj Seurat object (before filtering)
#' @param qc_params List of QC thresholds
#' @return Data frame summarizing filter failures
analyze_qc_failures <- function(obj, qc_params) {
  meta <- obj@meta.data
  
  failures <- tibble(
    cell = rownames(meta),
    fail_mito = meta$percent.mito > qc_params$mt_pct_max,
    fail_nfeature_low = meta$nFeature_RNA < qc_params$nfeature_min,
    fail_nfeature_high = meta$nFeature_RNA > qc_params$nfeature_max,
    fail_ncount_low = meta$nCount_RNA < qc_params$ncount_min,
    fail_ncount_high = meta$nCount_RNA > qc_params$ncount_max
  ) %>%
    mutate(
      n_failures = fail_mito + fail_nfeature_low + fail_nfeature_high + 
        fail_ncount_low + fail_ncount_high,
      pass = n_failures == 0
    )
  
  summary_stats <- failures %>%
    summarise(
      total_cells = n(),
      passed = sum(pass),
      failed = sum(!pass),
      pct_pass = 100 * passed / total_cells,
      fail_mito = sum(fail_mito),
      fail_low_features = sum(fail_nfeature_low),
      fail_high_features = sum(fail_nfeature_high),
      fail_low_counts = sum(fail_ncount_low),
      fail_high_counts = sum(fail_ncount_high)
    )
  
  return(list(per_cell = failures, summary = summary_stats))
}

# ============================================================================
# 3. GENE TRAJECTORY FUNCTIONS
# ============================================================================

#' Plot gene expression trajectories for paired samples
#' @param obj Seurat object
#' @param genes Vector of gene names
#' @param subject_col Subject ID column
#' @param time_col Treatment/time column
#' @param cluster_col Optional: cell type to subset
#' @return ggplot object
plot_gene_trajectory <- function(obj, genes, subject_col, time_col, 
                                 cluster_col = NULL, cluster_value = NULL) {
  
  # Subset to specific cluster if requested
  if (!is.null(cluster_col) && !is.null(cluster_value)) {
    cells <- rownames(obj@meta.data[obj@meta.data[[cluster_col]] == cluster_value, ])
    obj_sub <- subset(obj, cells = cells)
  } else {
    obj_sub <- obj
  }
  
  # Get expression data
  expr_data <- FetchData(obj_sub, vars = c(genes, subject_col, time_col))
  
  # Reshape to long format
  expr_long <- expr_data %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
    group_by(gene, !!sym(subject_col), !!sym(time_col)) %>%
    summarise(mean_expr = mean(expression), .groups = "drop")
  
  # Calculate paired statistics
  paired_stats <- expr_long %>%
    pivot_wider(names_from = time_col, values_from = mean_expr) %>%
    group_by(gene) %>%
    summarise(
      n_pairs = n(),
      mean_change = mean(get(names(.)[3]) - get(names(.)[4]), na.rm = TRUE),
      p_value = if (n() >= 3) {
        tryCatch(t.test(get(names(.)[3]), get(names(.)[4]), paired = TRUE)$p.value,
                 error = function(e) NA)
      } else {
        NA
      }
    ) %>%
    mutate(
      direction = ifelse(mean_change > 0, "Up", "Down"),
      sig_label = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  # Create plot
  p <- ggplot(expr_long, aes(x = !!sym(time_col), y = mean_expr, 
                             group = !!sym(subject_col))) +
    geom_line(aes(color = !!sym(subject_col)), alpha = 0.5, linewidth = 0.8) +
    geom_point(aes(color = !!sym(subject_col)), size = 2) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", 
                 color = "black", linewidth = 1.5, linetype = "dashed") +
    stat_summary(aes(group = 1), fun = mean, geom = "point", 
                 color = "black", size = 4, shape = 18) +
    facet_wrap(~ gene, scales = "free_y", ncol = min(4, length(genes))) +
    geom_text(data = paired_stats, 
              aes(x = 1.5, y = Inf, label = sig_label),
              vjust = 1.5, size = 6, inherit.aes = FALSE) +
    labs(
      title = "Gene Expression Trajectories (Paired Samples)",
      subtitle = if (!is.null(cluster_value)) paste("Cell type:", cluster_value) else "All cells",
      x = "Condition",
      y = "Mean Expression",
      color = "Subject"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold")
    )
  
  return(list(plot = p, statistics = paired_stats))
}

# ============================================================================
# 4. ENHANCED VOLCANO PLOTS
# ============================================================================

#' Create enhanced volcano plot with functional annotations
#' @param deg_results Data frame with DEG results
#' @param fdr_cutoff FDR threshold
#' @param fc_cutoff Log2FC threshold
#' @param top_n Number of top genes to label
#' @param gene_categories Optional named list of gene categories
#' @return ggplot object
plot_enhanced_volcano <- function(deg_results, fdr_cutoff = 0.05, fc_cutoff = 0.25,
                                  top_n = 20, gene_categories = NULL) {
  
  # Prepare data
  plot_data <- deg_results %>%
    mutate(
      log10_padj = -log10(adj.P.Val),
      significant = adj.P.Val <= fdr_cutoff & abs(logFC) >= fc_cutoff,
      direction = case_when(
        significant & logFC > 0 ~ "Up",
        significant & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Add gene categories if provided
  if (!is.null(gene_categories)) {
    for (cat_name in names(gene_categories)) {
      plot_data <- plot_data %>%
        mutate(category = ifelse(gene %in% gene_categories[[cat_name]], 
                                 cat_name, category))
    }
  }
  
  # Select top genes to label
  top_genes <- plot_data %>%
    filter(significant) %>%
    arrange(desc(abs(logFC))) %>%
    head(top_n)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = logFC, y = log10_padj)) +
    geom_point(aes(color = direction), alpha = 0.5, size = 1.5) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), 
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff), 
               linetype = "dashed", color = "grey50") +
    scale_color_manual(
      values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey70"),
      name = "Regulation"
    ) +
    labs(
      title = "Volcano Plot: Differential Expression",
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"Adjusted P-value")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  # Add gene labels
  if (nrow(top_genes) > 0) {
    p <- p + 
      geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5,
        segment.alpha = 0.5
      )
  }
  
  # Add category legend if provided
  if (!is.null(gene_categories)) {
    p <- p + 
      geom_point(data = plot_data %>% filter(!is.na(category)),
                 aes(shape = category), size = 3) +
      scale_shape_manual(values = c(15:18), name = "Gene Category")
  }
  
  return(p)
}

#' Create interactive volcano plot
#' @param deg_results Data frame with DEG results
#' @param fdr_cutoff FDR threshold
#' @param fc_cutoff Log2FC threshold
#' @return plotly object
plot_volcano_interactive <- function(deg_results, fdr_cutoff = 0.05, fc_cutoff = 0.25) {
  
  plot_data <- deg_results %>%
    mutate(
      log10_padj = -log10(adj.P.Val),
      significant = adj.P.Val <= fdr_cutoff & abs(logFC) >= fc_cutoff,
      direction = case_when(
        significant & logFC > 0 ~ "Up-regulated",
        significant & logFC < 0 ~ "Down-regulated",
        TRUE ~ "Not significant"
      ),
      hover_text = paste0(
        "Gene: ", gene, "<br>",
        "log2FC: ", round(logFC, 3), "<br>",
        "Adj. P-value: ", formatC(adj.P.Val, format = "e", digits = 2), "<br>",
        "Direction: ", direction
      )
    )
  
  p <- plot_ly(data = plot_data,
               x = ~logFC,
               y = ~log10_padj,
               color = ~direction,
               colors = c("Down-regulated" = "#4DBBD5", 
                          "Not significant" = "grey70",
                          "Up-regulated" = "#E64B35"),
               type = "scatter",
               mode = "markers",
               marker = list(size = 4, opacity = 0.6),
               text = ~hover_text,
               hoverinfo = "text") %>%
    add_segments(x = -fc_cutoff, xend = -fc_cutoff, y = 0, yend = max(plot_data$log10_padj),
                 line = list(dash = "dash", color = "grey"), 
                 showlegend = FALSE, hoverinfo = "none") %>%
    add_segments(x = fc_cutoff, xend = fc_cutoff, y = 0, yend = max(plot_data$log10_padj),
                 line = list(dash = "dash", color = "grey"), 
                 showlegend = FALSE, hoverinfo = "none") %>%
    add_segments(x = min(plot_data$logFC), xend = max(plot_data$logFC),
                 y = -log10(fdr_cutoff), yend = -log10(fdr_cutoff),
                 line = list(dash = "dash", color = "grey"), 
                 showlegend = FALSE, hoverinfo = "none") %>%
    layout(
      title = "Interactive Volcano Plot",
      xaxis = list(title = "log2 Fold Change"),
      yaxis = list(title = "-log10 Adjusted P-value"),
      hovermode = "closest"
    )
  
  return(p)
}

# ============================================================================
# 5. LEADING EDGE ANALYSIS
# ============================================================================

#' Extract and visualize leading edge genes from GSEA
#' @param gsea_result GSEA result object
#' @param deg_results DEG results data frame
#' @param pathway_id Pathway ID to visualize
#' @param top_n Number of top pathways to show
#' @return List with plot and leading edge genes
plot_leading_edge <- function(gsea_result, deg_results, pathway_id = NULL, top_n = 10) {
  
  # If no pathway specified, take top pathways
  if (is.null(pathway_id)) {
    top_pathways <- gsea_result %>%
      arrange(p.adjust) %>%
      head(top_n) %>%
      pull(ID)
  } else {
    top_pathways <- pathway_id
  }
  
  # Extract leading edge genes
  leading_edge_list <- list()
  for (pw in top_pathways) {
    pw_data <- gsea_result %>% filter(ID == pw)
    if (nrow(pw_data) > 0 && !is.null(pw_data$core_enrichment)) {
      genes <- strsplit(pw_data$core_enrichment, "/")[[1]]
      leading_edge_list[[pw]] <- genes
    }
  }
  
  # Create heatmap data
  all_le_genes <- unique(unlist(leading_edge_list))
  
  expr_matrix <- deg_results %>%
    filter(gene %in% all_le_genes) %>%
    select(gene, logFC) %>%
    column_to_rownames("gene")
  
  # Add pathway membership annotations
  pathway_membership <- matrix(0, nrow = length(all_le_genes), ncol = length(leading_edge_list))
  rownames(pathway_membership) <- all_le_genes
  colnames(pathway_membership) <- names(leading_edge_list)
  
  for (i in seq_along(leading_edge_list)) {
    pw_name <- names(leading_edge_list)[i]
    genes <- leading_edge_list[[i]]
    pathway_membership[genes, pw_name] <- 1
  }
  
  return(list(
    leading_edge_genes = leading_edge_list,
    expression_matrix = expr_matrix,
    pathway_membership = pathway_membership
  ))
}

# ============================================================================
# 6. GENE SIGNATURE SCORING
# ============================================================================

#' Score cells for custom gene signatures
#' @param obj Seurat object
#' @param signatures Named list of gene vectors
#' @param method Scoring method: "Seurat" or "UCell"
#' @return Seurat object with scores added to metadata
score_gene_signatures <- function(obj, signatures, method = "Seurat") {
  
  log_info("Scoring {length(signatures)} gene signatures using {method} method")
  
  if (method == "Seurat") {
    for (sig_name in names(signatures)) {
      genes <- intersect(signatures[[sig_name]], rownames(obj))
      if (length(genes) < 3) {
        log_warn("Signature '{sig_name}' has fewer than 3 genes in dataset, skipping")
        next
      }
      
      obj <- AddModuleScore(
        obj,
        features = list(genes),
        name = paste0(sig_name, "_score"),
        assay = "RNA"
      )
      
      # Rename column (Seurat adds "1" to name)
      colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0(sig_name, "_score")
    }
  } else if (method == "UCell") {
    if (requireNamespace("UCell", quietly = TRUE)) {
      obj <- UCell::AddModuleScore_UCell(obj, features = signatures)
    } else {
      log_warn("UCell not installed, falling back to Seurat method")
      return(score_gene_signatures(obj, signatures, method = "Seurat"))
    }
  }
  
  log_success("Gene signature scoring complete")
  return(obj)
}

# ============================================================================
# 7. SAMPLE-LEVEL SUMMARY
# ============================================================================

#' Generate comprehensive sample quality report
#' @param obj Seurat object
#' @param subject_col Subject ID column
#' @param treatment_col Treatment column
#' @return Data frame with per-sample metrics
generate_sample_report <- function(obj, subject_col, treatment_col) {
  
  sample_stats <- obj@meta.data %>%
    group_by(!!sym(subject_col), !!sym(treatment_col)) %>%
    summarise(
      n_cells = n(),
      median_UMI = median(nCount_RNA),
      median_genes = median(nFeature_RNA),
      median_mito = median(percent.mito),
      median_ribo = median(percent.ribo, na.rm = TRUE),
      mean_novelty = mean(nFeature_RNA / nCount_RNA, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(!!sym(subject_col), !!sym(treatment_col))
  
  # Flag outliers
  sample_stats <- sample_stats %>%
    mutate(
      flag_low_cells = n_cells < quantile(n_cells, 0.25) - 1.5 * IQR(n_cells),
      flag_high_mito = median_mito > quantile(median_mito, 0.75) + 1.5 * IQR(median_mito),
      flag_low_genes = median_genes < quantile(median_genes, 0.25) - 1.5 * IQR(median_genes),
      outlier = flag_low_cells | flag_high_mito | flag_low_genes
    )
  
  return(sample_stats)
}

# ============================================================================
# 8. EXPORT FUNCTIONS
# ============================================================================

#' Export data for external visualization tools
#' @param obj Seurat object
#' @param output_dir Output directory
#' @param formats Vector of formats: "cellxgene", "loupe", "ucsc"
export_for_viewers <- function(obj, output_dir, formats = c("cellxgene")) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if ("cellxgene" %in% formats) {
    log_info("Exporting for cellxgene...")
    # Save as h5ad format
    if (requireNamespace("SeuratDisk", quietly = TRUE)) {
      h5seurat_file <- file.path(output_dir, "data.h5Seurat")
      h5ad_file <- file.path(output_dir, "data.h5ad")
      
      SeuratDisk::SaveH5Seurat(obj, filename = h5seurat_file)
      SeuratDisk::Convert(h5seurat_file, dest = "h5ad")
      
      log_success("cellxgene export complete: {h5ad_file}")
    } else {
      log_warn("SeuratDisk not installed, skipping cellxgene export")
    }
  }
  
  if ("loupe" %in% formats) {
    log_info("Exporting for Loupe Browser...")
    # Save as Loupe file (requires loupeR package)
    if (requireNamespace("loupeR", quietly = TRUE)) {
      loupe_file <- file.path(output_dir, "data.cloupe")
      loupeR::create_loupe_from_seurat(obj, output_name = loupe_file)
      log_success("Loupe export complete: {loupe_file}")
    } else {
      log_warn("loupeR not installed, skipping Loupe export")
    }
  }
  
  invisible(NULL)
}

#' Add hyperlinks to gene names for external databases
#' @param gene_names Vector of gene names
#' @param database Database to link to: "genecards", "ensembl", "ncbi"
#' @return Vector of HTML links
add_gene_links <- function(gene_names, database = "genecards") {
  
  urls <- switch(database,
                 "genecards" = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene_names),
                 "ensembl" = paste0("https://www.ensembl.org/Multi/Search/Results?q=", gene_names),
                 "ncbi" = paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", gene_names),
                 stop("Unknown database")
  )
  
  html_links <- sprintf('<a href="%s" target="_blank">%s</a>', urls, gene_names)
  return(html_links)
}

# ============================================================================
# 9. PROGRESS TRACKING
# ============================================================================

#' Create progress bar for long operations
#' @param total Total number of iterations
#' @param format Progress bar format
#' @return progress_bar object
create_progress_bar <- function(total, format = "[:bar] :percent ETA: :eta") {
  pb <- progress_bar$new(
    format = format,
    total = total,
    clear = FALSE,
    width = 60
  )
  return(pb)
}

# ============================================================================
# 10. UTILITY FUNCTIONS
# ============================================================================

#' Save plot in multiple formats
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param output_dir Output directory
#' @param formats Vector of formats: "png", "pdf", "svg"
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi DPI for raster formats
save_plot_multi <- function(plot, filename, output_dir, 
                            formats = c("png", "pdf"), 
                            width = 10, height = 8, dpi = 300) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (fmt in formats) {
    file_path <- file.path(output_dir, paste0(filename, ".", fmt))
    
    ggsave(
      filename = file_path,
      plot = plot,
      width = width,
      height = height,
      dpi = if (fmt %in% c("png", "jpg")) dpi else NULL,
      device = fmt
    )
    
    log_info("Saved: {file_path}")
  }
  
  invisible(NULL)
}

#' Format p-values for display
#' @param p_values Numeric vector of p-values
#' @return Character vector of formatted p-values
format_pvalue <- function(p_values) {
  sapply(p_values, function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("< 0.001")
    if (p < 0.01) return(sprintf("%.3f", p))
    return(sprintf("%.2f", p))
  })
}

log_info("Enhanced helper functions loaded successfully")