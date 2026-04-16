#PBMC_DEG_ROOT  <- "deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned/paired_cross_pool/deg_by_celltype"   # folder of per-celltype DEG CSVs
#ILEUM_DEG_ROOT <- "/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/OUTPUT_20251128_141815_WITH_BATCH/Monaco_fine/deg_results"  # folder of per-celltype DEG CSVs

#PBMC_DEG_ROOT  <- "deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned/paired_cross_pool/deg_by_celltype"
#ILEUM_DEG_ROOT <- "/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/OUTPUT_20251128_141815_WITH_BATCH/Monaco_fine/deg_results"

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(fs)
  library(purrr)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

SCATTER_GENES_TO_SHOW <- c(
  # interferon / antiviral
  "STAT1", "IRF1", "IRF9",
  "ISG15", "CXCL10",
  "GBP1", "GBP5",
  "PARP9", "PLSCR1",
  "TAP1", "TAP2",
  "PSMB8", "PSMB9",
  "NLRC5",
  
  # antigen presentation
  "HLA-DPA1", "HLA-DPB1", "HLA-DRB1", "CD74",
  
  # cytotoxic / T cell
  "GZMB", "GZMA", "GZMH",
  "IL7R", "CXCR6", "KLRB1",
  
  # B cell / plasma
  "JCHAIN", "MZB1", "TNFRSF17",
  
  # signaling / survival
  "BCL2", "AKT1", "VAV3"
)

PBMC_DEG_ROOT  <- "deg_gsea_pooled_vs_within_pool_Monaco_fine_pruned/paired_cross_pool/deg_by_celltype"
ILEUM_DEG_ROOT <- "/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/OUTPUT_20251128_141815_WITH_BATCH/Monaco_fine/deg_results"

OUT_DIR <- "scatter_gene_outputs"
dir_create(OUT_DIR)

FDR_STRICT    <- 0.05
FDR_LOOSE     <- 0.10
TOP_N_SCATTER <- 3
LABEL_WIDTH   <- 20
MIN_ABS_LFC   <- 0

PBMC_COL_MAP <- c(
  gene           = "gene",
  log2FoldChange = "logFC",
  padj           = "adj.P.Val"
)

ILEUM_COL_MAP <- c(
  gene           = "gene",
  log2FoldChange = "logFC",
  padj           = "adj.P.Val"
)

# =============================================================================
# HELPERS
# =============================================================================

clean_cluster_names <- function(x) {
  x %>%
    str_replace_all("\\+", "") %>%
    str_replace_all("__", "_") %>%
    str_replace_all(" ", "_") %>%
    str_replace_all("-", "_") %>%
    str_replace_all("\\.", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace("^cluster_", "") %>%
    str_replace("^DEG_", "") %>%
    str_replace("_DEG$", "") %>%
    str_replace("^ALL_clusters_", "") %>%
    str_replace("^_|_$", "") %>%
    str_trim()
}

save_plot <- function(plot_obj, filename, width, height, dpi = 300) {
  out <- file.path(OUT_DIR, filename)
  ggsave(out, plot_obj, width = width, height = height, dpi = dpi, limitsize = FALSE)
  invisible(out)
}

remap_cols <- function(df, col_map) {
  for (std_name in names(col_map)) {
    raw_name <- col_map[[std_name]]
    if (raw_name %in% colnames(df) && raw_name != std_name) {
      df <- df %>% rename(!!std_name := !!raw_name)
    }
  }
  df
}

extract_celltype_from_path <- function(f, root) {
  x <- basename(f) %>% str_remove("\\.csv$")
  
  x %>%
    str_remove("^Monaco_fine_cluster_") %>%
    str_remove("^Monaco_fine_") %>%
    str_remove("^DEG_") %>%
    str_remove("_DEG$") %>%
    str_remove("_limma_Treated_vs_Untreated$") %>%
    str_remove("_dream_Treated_vs_Untreated$") %>%
    str_remove("_Treated_vs_Untreated$") %>%
    str_remove("_MildTreated_vs_MildUntreated.*$") %>%
    str_remove("_results.*$") %>%
    clean_cluster_names()
}

dedup_deg <- function(df) {
  df %>%
    group_by(gene, cluster, dataset) %>%
    summarise(
      log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))][1],
      padj = {
        x <- padj[!is.na(padj)]
        if (length(x) == 0) NA_real_ else min(x)
      },
      .groups = "drop"
    )
}

# =============================================================================
# 1. LOAD PBMC
# =============================================================================

pbmc_files <- dir_ls(PBMC_DEG_ROOT, recurse = TRUE, regexp = "\\.csv$")

load_pbmc_deg <- function(f) {
  df <- fread(f, fill = TRUE)
  
  if (!"gene" %in% colnames(df) && !"logFC" %in% colnames(df)) {
    message("Skipping PBMC file with no DEG columns: ", basename(f))
    return(NULL)
  }
  
  df <- remap_cols(df, PBMC_COL_MAP)
  df$cluster <- extract_celltype_from_path(f, PBMC_DEG_ROOT)
  df$dataset <- "PBMC"
  df
}

pbmc_deg <- map(pbmc_files, load_pbmc_deg) %>%
  keep(~ !is.null(.)) %>%
  bind_rows()

# =============================================================================
# 2. LOAD ILEUM
# Keep only per-cluster limma files; exclude dream and ALL_clusters
# =============================================================================

ileum_files <- dir_ls(ILEUM_DEG_ROOT, recurse = TRUE, regexp = "\\.csv$") %>%
  keep(~ str_detect(basename(.x), "limma")) %>%
  keep(~ str_detect(basename(.x), "cluster")) %>%
  keep(~ !str_detect(basename(.x), "ALL_clusters"))

load_ileum_deg <- function(f) {
  df <- fread(f, fill = TRUE)
  
  if (!"gene" %in% colnames(df) && !"logFC" %in% colnames(df)) {
    message("Skipping ILEUM file with no DEG columns: ", basename(f))
    return(NULL)
  }
  
  df <- remap_cols(df, ILEUM_COL_MAP)
  df$cluster <- extract_celltype_from_path(f, ILEUM_DEG_ROOT)
  df$dataset <- "ILEUM"
  df
}

ileum_deg <- map(ileum_files, load_ileum_deg) %>%
  keep(~ !is.null(.)) %>%
  bind_rows()

# =============================================================================
# 3. CLEAN + COMBINE
# =============================================================================

pbmc_deg <- pbmc_deg %>%
  mutate(
    gene = toupper(trimws(gene)),
    cluster = clean_cluster_names(cluster),
    dataset = "PBMC"
  )

ileum_deg <- ileum_deg %>%
  mutate(
    gene = toupper(trimws(gene)),
    cluster = clean_cluster_names(cluster),
    dataset = "ILEUM"
  )

deg_all <- bind_rows(pbmc_deg, ileum_deg) %>%
  filter(!is.na(gene), !is.na(log2FoldChange), cluster != "")

if (MIN_ABS_LFC > 0) {
  deg_all <- deg_all %>% filter(abs(log2FoldChange) >= MIN_ABS_LFC)
}

cat("PBMC rows:", nrow(pbmc_deg), "\n")
cat("ILEUM rows:", nrow(ileum_deg), "\n")
cat("Datasets in deg_all:\n")
print(deg_all %>% count(dataset))

cat("Shared clusters:\n")
print(intersect(
  sort(unique(pbmc_deg$cluster)),
  sort(unique(ileum_deg$cluster))
))

# =============================================================================
# 4. DEDUP + WIDEN
# =============================================================================

datasets_present <- unique(deg_all$dataset)

if (!"PBMC" %in% datasets_present) {
  stop("PBMC rows are missing from deg_all")
}
if (!"ILEUM" %in% datasets_present) {
  stop("ILEUM rows are missing from deg_all. Check ileum_files filtering.")
}

wide_deg <- deg_all %>%
  select(gene, cluster, dataset, log2FoldChange, padj) %>%
  dedup_deg() %>%
  pivot_wider(
    names_from  = dataset,
    values_from = c(log2FoldChange, padj),
    names_sep   = "__"
  )

required_cols <- c(
  "log2FoldChange__PBMC",
  "log2FoldChange__ILEUM",
  "padj__PBMC",
  "padj__ILEUM"
)

missing_cols <- setdiff(required_cols, colnames(wide_deg))
if (length(missing_cols) > 0) {
  stop("Missing expected columns after pivot_wider(): ", paste(missing_cols, collapse = ", "))
}

matched_genes <- wide_deg %>%
  filter(!is.na(log2FoldChange__PBMC), !is.na(log2FoldChange__ILEUM)) %>%
  mutate(
    sig_status = case_when(
      !is.na(padj__PBMC) & !is.na(padj__ILEUM) &
        padj__PBMC <= FDR_LOOSE & padj__ILEUM <= FDR_LOOSE ~ "Both significant",
      !is.na(padj__PBMC) & padj__PBMC <= FDR_LOOSE &
        (is.na(padj__ILEUM) | padj__ILEUM > FDR_LOOSE) ~ "PBMC only",
      !is.na(padj__ILEUM) & padj__ILEUM <= FDR_LOOSE &
        (is.na(padj__PBMC) | padj__PBMC > FDR_LOOSE) ~ "ILEUM only",
      TRUE ~ "Neither significant"
    ),
    sig_tier = case_when(
      !is.na(padj__PBMC) & !is.na(padj__ILEUM) &
        padj__PBMC <= FDR_STRICT & padj__ILEUM <= FDR_STRICT ~ "Both ??? 0.05",
      
      !is.na(padj__PBMC) & !is.na(padj__ILEUM) &
        padj__PBMC <= FDR_LOOSE & padj__ILEUM <= FDR_LOOSE ~ "Both / mixed 0.05???0.1",
      
      !is.na(padj__PBMC) & padj__PBMC <= FDR_STRICT &
        (is.na(padj__ILEUM) | padj__ILEUM > FDR_LOOSE) ~ "PBMC only ??? 0.05",
      
      !is.na(padj__PBMC) & padj__PBMC > FDR_STRICT & padj__PBMC <= FDR_LOOSE &
        (is.na(padj__ILEUM) | padj__ILEUM > FDR_LOOSE) ~ "PBMC only 0.05???0.1",
      
      !is.na(padj__ILEUM) & padj__ILEUM <= FDR_STRICT &
        (is.na(padj__PBMC) | padj__PBMC > FDR_LOOSE) ~ "ILEUM only ??? 0.05",
      
      !is.na(padj__ILEUM) & padj__ILEUM > FDR_STRICT & padj__ILEUM <= FDR_LOOSE &
        (is.na(padj__PBMC) | padj__PBMC > FDR_LOOSE) ~ "ILEUM only 0.05???0.1",
      
      TRUE ~ "Not significant"
    ),
    quadrant = case_when(
      log2FoldChange__PBMC >  0 & log2FoldChange__ILEUM >  0 ~ "Q1: Both Up",
      log2FoldChange__PBMC <  0 & log2FoldChange__ILEUM >  0 ~ "Q2: Ileum Up / PBMC Down",
      log2FoldChange__PBMC <  0 & log2FoldChange__ILEUM <  0 ~ "Q3: Both Down",
      log2FoldChange__PBMC >  0 & log2FoldChange__ILEUM <  0 ~ "Q4: PBMC Up / Ileum Down",
      TRUE ~ "On Axis"
    ),
    label_strength = abs(log2FoldChange__PBMC) + abs(log2FoldChange__ILEUM)
  )

fwrite(matched_genes, file.path(OUT_DIR, "matched_genes_PBMC_vs_ILEUM.csv"))

cat("Matched genes:", nrow(matched_genes), "\n")
cat("Matched clusters:", length(unique(matched_genes$cluster)), "\n")

# =============================================================================
# 5. LABEL SETUP
# =============================================================================

scatter_df <- matched_genes %>%
  filter(sig_status != "Neither significant", !is.na(cluster))

top_genes_auto <- if (nrow(scatter_df) > 0) {
  scatter_df %>%
    group_by(cluster, quadrant) %>%
    slice_max(order_by = label_strength, n = TOP_N_SCATTER, with_ties = FALSE) %>%
    ungroup()
} else {
  scatter_df %>% slice(0)
}

top_genes_user <- if (!is.null(SCATTER_GENES_TO_SHOW) && length(SCATTER_GENES_TO_SHOW) > 0) {
  scatter_df %>% filter(gene %in% toupper(SCATTER_GENES_TO_SHOW))
} else {
  scatter_df %>% slice(0)
}

top_genes <- bind_rows(top_genes_auto, top_genes_user) %>%
  distinct(cluster, gene, .keep_all = TRUE)

# =============================================================================
# 6. PLOTS
# =============================================================================

quadrant_colors <- c(
  "Q1: Both Up"              = "#D6604D",
  "Q2: Ileum Up / PBMC Down" = "#4393C3",
  "Q3: Both Down"            = "#2166AC",
  "Q4: PBMC Up / Ileum Down" = "#B2182B",
  "On Axis"                  = "grey70"
)

sig_shapes <- c(
  "Both significant" = 16,
  "PBMC only"        = 17,
  "ILEUM only"       = 15
)

sig_alpha <- c(
  "Both ??? 0.05"            = 0.95,
  "Both / mixed 0.05???0.1"  = 0.45,
  "PBMC only ??? 0.05"       = 0.95,
  "PBMC only 0.05???0.1"     = 0.45,
  "ILEUM only ??? 0.05"      = 0.95,
  "ILEUM only 0.05???0.1"    = 0.45
)

if (nrow(scatter_df) == 0) {
  message("No significant matched genes to plot at current FDR cutoff.")
} else {
  p_scatter_genes <- ggplot(
    scatter_df,
    aes(
      x = log2FoldChange__PBMC,
      y = log2FoldChange__ILEUM,
      color = quadrant,
      shape = sig_status,
      alpha = sig_tier
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_point(size = 1.1) +
    geom_text_repel(
      data = top_genes,
      aes(label = str_trunc(gene, LABEL_WIDTH)),
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
    scale_alpha_manual(
      values = sig_alpha,
      breaks = c(
        "Both ??? 0.05",
        "Both / mixed 0.05???0.1",
        "PBMC only ??? 0.05",
        "PBMC only 0.05???0.1",
        "ILEUM only ??? 0.05",
        "ILEUM only 0.05???0.1"
      ),
      name = "FDR tier"
    ) +
    labs(
      title = "PBMC vs Ileum gene log2FC correlation",
      subtitle = paste0(
        "Genes significant at FDR ??? ", FDR_LOOSE,
        "; darker = ??? ", FDR_STRICT,
        ", lighter = 0.05???0.1; top ",
        TOP_N_SCATTER, " per quadrant labelled"
      ),
      x = "PBMC log2FoldChange",
      y = "Ileum log2FoldChange",
      color = "Quadrant",
      shape = "Significance"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  save_plot(p_scatter_genes, "scatter_genes_PBMC_vs_ILEUM.pdf", 20, 16)
  save_plot(p_scatter_genes, "scatter_genes_PBMC_vs_ILEUM.png", 20, 16)
}

discord_df <- matched_genes %>%
  filter(
    sig_status != "Neither significant",
    quadrant %in% c("Q2: Ileum Up / PBMC Down", "Q4: PBMC Up / Ileum Down")
  )

if (nrow(discord_df) == 0) {
  message("No discordant significant genes to plot.")
} else {
  p_scatter_genes_discord <- ggplot(
    discord_df,
    aes(
      x = log2FoldChange__PBMC,
      y = log2FoldChange__ILEUM,
      color = quadrant,
      shape = sig_status,
      alpha = sig_tier
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_point(size = 1.4) +
    facet_wrap(~ cluster, scales = "free", ncol = 3) +
    scale_color_manual(values = quadrant_colors) +
    scale_shape_manual(values = sig_shapes) +
    scale_alpha_manual(
      values = sig_alpha,
      breaks = c(
        "Both ??? 0.05",
        "Both / mixed 0.05???0.1",
        "PBMC only ??? 0.05",
        "PBMC only 0.05???0.1",
        "ILEUM only ??? 0.05",
        "ILEUM only 0.05???0.1"
      ),
      name = "FDR tier"
    ) +
    labs(
      title = "PBMC vs Ileum discordant genes",
      subtitle = paste0(
        "Genes significant at FDR ??? ", FDR_LOOSE,
        "; darker = ??? ", FDR_STRICT,
        ", lighter = 0.05???0.1"
      ),
      x = "PBMC log2FoldChange",
      y = "Ileum log2FoldChange",
      color = "Quadrant",
      shape = "Significance"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  save_plot(p_scatter_genes_discord, "scatter_genes_PBMC_vs_ILEUM_discordant_only.pdf", 18, 12)
}

# =============================================================================
# 7. QUADRANT SUMMARY
# =============================================================================

if (nrow(matched_genes) > 0) {
  gene_quad_counts <- matched_genes %>%
    count(cluster, quadrant, sort = TRUE)
  
  fwrite(gene_quad_counts, file.path(OUT_DIR, "gene_quadrant_counts.csv"))
  
  p_gene_quad_counts <- ggplot(gene_quad_counts, aes(x = cluster, y = n, fill = quadrant)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = quadrant_colors) +
    labs(
      title = "Number of matched genes by scatter quadrant",
      x = "Cell type",
      y = "Number of genes"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_plot(p_gene_quad_counts, "summary_gene_quadrant_counts.pdf", 11, 6)
}

cat("Done. Outputs written to:", OUT_DIR, "\n")
cat("Matched genes (PBMC + ILEUM):", nrow(matched_genes), "\n")
cat("Significant in at least one tissue:", nrow(scatter_df), "\n")

sort(unique(pbmc_deg$cluster))
sort(unique(ileum_deg$cluster))
intersect(sort(unique(pbmc_deg$cluster)), sort(unique(ileum_deg$cluster)))