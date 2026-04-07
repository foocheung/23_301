# ── 1. Per-subject DEG consistency check ─────────────────────────────────────
# For each subject, run a simple condition comparison within monocytes only
# Then ask: are the same genes going in the same direction across subjects?

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)

obj <- readRDS("PBMC_RDS/object_final_annotated.rds")

SUBJECT_COL    <- "Subject_raw"
CONDITION_COL  <- "Treatment"
CLUSTER_COL_L1 <- "Monaco_main"
CONTRAST       <- c("Treated", "Untreated")   # adjust to your actual labels

# Subset to monocytes only
mono_cells <- rownames(obj@meta.data)[obj@meta.data[[CLUSTER_COL_L1]] == "Monocytes"]
obj_mono   <- subset(obj, cells = mono_cells)

subjects_with_both <- obj_mono@meta.data %>%
  group_by(.data[[SUBJECT_COL]]) %>%
  summarise(
    conditions = list(unique(.data[[CONDITION_COL]])),
    .groups = "drop"
  ) %>%
  filter(sapply(conditions, function(x) all(CONTRAST %in% x))) %>%
  pull(.data[[SUBJECT_COL]])

message("Subjects with both conditions in monocytes: ",
        paste(subjects_with_both, collapse = ", "))

# Per-subject pseudo-bulk or simple FindMarkers
per_subject_degs <- lapply(subjects_with_both, function(subj) {
  cells_subj <- rownames(obj_mono@meta.data)[
    obj_mono@meta.data[[SUBJECT_COL]] == subj
  ]
  obj_subj <- subset(obj_mono, cells = cells_subj)
  
  # Need minimum cells per condition
  n_per_cond <- table(obj_subj@meta.data[[CONDITION_COL]])
  if (any(n_per_cond < 10)) {
    message(sprintf("Skipping %s — too few cells: %s", subj,
                    paste(names(n_per_cond), n_per_cond, sep = "=", collapse = ", ")))
    return(NULL)
  }
  
  Idents(obj_subj) <- CONDITION_COL
  tryCatch({
    markers <- FindMarkers(
      obj_subj,
      ident.1     = CONTRAST[1],
      ident.2     = CONTRAST[2],
      test.use    = "wilcox",
      min.pct     = 0.1,
      logfc.threshold = 0.1,
      verbose     = FALSE
    )
    markers$gene    <- rownames(markers)
    markers$subject <- subj
    markers
  }, error = function(e) {
    message(sprintf("Error for %s: %s", subj, e$message))
    NULL
  })
})

names(per_subject_degs) <- subjects_with_both
per_subject_degs        <- Filter(Negate(is.null), per_subject_degs)

# ── 2. LogFC concordance heatmap across subjects ──────────────────────────────
# Get top genes from your overall DEG result and check direction per subject
top_genes <- bind_rows(per_subject_degs) %>%
  group_by(gene) %>%
  summarise(
    n_subjects   = n(),
    mean_abs_lfc = mean(abs(avg_log2FC)),
    .groups = "drop"
  ) %>%
  filter(n_subjects >= max(1, length(per_subject_degs) - 1)) %>%  # in most subjects
  slice_max(mean_abs_lfc, n = 50) %>%
  pull(gene)

lfc_matrix <- bind_rows(per_subject_degs) %>%
  filter(gene %in% top_genes) %>%
  select(gene, subject, avg_log2FC) %>%
  tidyr::pivot_wider(names_from = subject, values_from = avg_log2FC, values_fill = 0) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Colour scale centred at 0
lfc_max <- max(abs(lfc_matrix), na.rm = TRUE)
breaks  <- seq(-lfc_max, lfc_max, length.out = 101)
colours <- colorRampPalette(c("#3C5488", "white", "#E64B35"))(100)

png("monocyte_persubject_logFC_heatmap.png",
    width = 10, height = 14, units = "in", res = 300, bg = "white")
pheatmap(
  lfc_matrix,
  color            = colours,
  breaks           = breaks,
  cluster_rows     = TRUE,
  cluster_cols     = FALSE,
  fontsize_row     = 7,
  fontsize_col     = 10,
  main             = sprintf("Monocyte DEG log2FC per subject\n%s vs %s",
                             CONTRAST[1], CONTRAST[2]),
  border_color     = NA
)
dev.off()
message("Saved: monocyte_persubject_logFC_heatmap.png")

# ── 3. Concordance summary — what fraction of subjects agree on direction? ────
concordance <- bind_rows(per_subject_degs) %>%
  filter(gene %in% top_genes) %>%
  group_by(gene) %>%
  summarise(
    n_up        = sum(avg_log2FC > 0),
    n_down      = sum(avg_log2FC < 0),
    n_subjects  = n(),
    pct_concordant = max(n_up, n_down) / n() * 100,
    mean_lfc    = mean(avg_log2FC),
    .groups     = "drop"
  ) %>%
  arrange(desc(pct_concordant))

message("\nTop concordant monocyte DEGs across subjects:")
print(head(concordance, 20))

png("monocyte_concordance_plot.png",
    width = 10, height = 8, units = "in", res = 300, bg = "white")
concordance %>%
  slice_max(mean_abs_lfc <- abs(mean_lfc), n = 40) %>%
  mutate(gene = reorder(gene, pct_concordant)) %>%
  ggplot(aes(x = pct_concordant, y = gene, fill = mean_lfc)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 75, linetype = "dashed", colour = "grey50") +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Mean log2FC") +
  labs(
    title    = "Monocyte DEG cross-subject concordance",
    subtitle = sprintf("%s vs %s · dashed line = 75%% concordance threshold",
                       CONTRAST[1], CONTRAST[2]),
    x        = "% subjects with concordant direction",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
dev.off()
message("Saved: monocyte_concordance_plot.png")

# ── 4. Variance partitioning — condition vs subject contribution ──────────────
# How much variance in key monocyte genes is explained by condition vs subject?
if (requireNamespace("variancePartition", quietly = TRUE)) {
  library(variancePartition)
  
  counts_mono <- GetAssayData(obj_mono, layer = "counts")
  meta_mono   <- obj_mono@meta.data
  
  # Keep only genes expressed in >10% of cells
  keep_genes <- rowSums(counts_mono > 0) / ncol(counts_mono) > 0.10
  counts_filt <- counts_mono[keep_genes, ]
  
  # Variance partition formula: condition + subject as random effects
  form <- as.formula(sprintf("~ (1|%s) + (1|%s)", CONDITION_COL, SUBJECT_COL))
  
  vp <- fitExtractVarPartModel(
    counts_filt,
    form,
    meta_mono,
    BPPARAM = BiocParallel::MulticoreParam(4)   # adjust cores as needed
  )
  
  png("monocyte_variance_partition.png",
      width = 10, height = 6, units = "in", res = 300, bg = "white")
  plotVarPart(sortCols(vp),
              main = "Variance partition — Monocytes\nCondition vs Subject contribution")
  dev.off()
  message("Saved: monocyte_variance_partition.png")
  
  # Print summary for interferon / JAK-STAT genes specifically
  ifn_genes <- intersect(
    c("IFIT1","IFIT2","IFIT3","ISG15","MX1","OAS1","OAS2",
      "STAT1","STAT2","IRF7","IRF9","IFI44L","IFI6","RSAD2"),
    rownames(vp)
  )
  if (length(ifn_genes) > 0) {
    message("\nVariance partition for IFN/JAK-STAT genes:")
    print(vp[ifn_genes, ])
  }
} else {
  message("variancePartition not installed — skipping variance decomposition")
  message("Install with: BiocManager::install('variancePartition')")
}

message("\n── All outputs ──────────────────────────────────────────────────────")
message("  monocyte_persubject_logFC_heatmap.png  — DEG direction per subject")
message("  monocyte_concordance_plot.png          — % subjects agreeing on direction")
message("  monocyte_variance_partition.png        — condition vs subject variance")