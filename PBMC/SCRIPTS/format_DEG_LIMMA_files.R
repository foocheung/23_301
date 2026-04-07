suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(fs)
})

# =========================
# USER VARIABLES
# =========================
RDS_FILE <- "./PBMC_RDS/object_final_annotated.rds"

CLUSTER_COL  <- "Monaco_fine"
CLUSTER_COL2 <- "Monaco_main"

SRC_BASE <- "deg_gsea_pooled_vs_within_pool_v2"
OUT_ROOT <- "./PBMC_report_outputs"

analysis_types <- c("paired_cross_pool", "pooled")
cluster_cols   <- c(CLUSTER_COL, CLUSTER_COL2)

# =========================
# HELPERS
# =========================
safe_fread <- function(f) {
  tryCatch(
    data.table::fread(f),
    error = function(e) {
      message("Skipping unreadable file: ", f)
      NULL
    }
  )
}

get_celltype <- function(f) basename(dirname(f))

bind_deg_files <- function(files, cluster_col_name) {
  if (length(files) == 0) return(NULL)
  
  purrr::map_dfr(files, function(f) {
    x <- safe_fread(f)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    
    x <- as.data.frame(x)
    
    if (!"cluster" %in% names(x)) {
      x$cluster <- get_celltype(f)
    }
    
    x$cluster_source <- cluster_col_name
    x$source_file <- basename(f)
    x$source_path <- f
    x
  })
}

bind_gsea_files <- function(files, cluster_col_name, geneset_name) {
  if (length(files) == 0) return(NULL)
  
  purrr::map_dfr(files, function(f) {
    x <- safe_fread(f)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    
    x <- as.data.frame(x)
    
    if (!"cluster" %in% names(x)) {
      x$cluster <- get_celltype(f)
    }
    
    x$cluster_source <- cluster_col_name
    x$geneset <- geneset_name
    x$source_file <- basename(f)
    x$source_path <- f
    x
  })
}

write_csv_if_present <- function(df, out) {
  if (is.null(df) || nrow(df) == 0) {
    message("No data found for: ", out)
    return(invisible(NULL))
  }
  dir_create(path_dir(out), recurse = TRUE)
  data.table::fwrite(df, out)
  message("Wrote: ", out, " [", nrow(df), " rows]")
}

make_combined_set <- function(analysis_type, cluster_col_name) {
  message("\n==============================")
  message("Building: ", analysis_type, " / ", cluster_col_name)
  message("==============================")
  
  src_deg  <- file.path(SRC_BASE, analysis_type, "deg_by_celltype")
  src_gsea <- file.path(SRC_BASE, analysis_type, "gsea_by_celltype")
  
  out_deg  <- file.path(OUT_ROOT, analysis_type, cluster_col_name, "deg_results")
  out_gsea <- file.path(OUT_ROOT, analysis_type, cluster_col_name, "gsea_results")
  
  limma_path <- file.path(out_deg,  sprintf("%s__ALL_clusters_limma.csv", cluster_col_name))
  go_path    <- file.path(out_gsea, sprintf("%s__ALL_clusters_GSEA_GO.csv", cluster_col_name))
  re_path    <- file.path(out_gsea, sprintf("%s__ALL_clusters_GSEA_REACTOME.csv", cluster_col_name))
  
  deg_files <- list.files(
    src_deg,
    pattern = "^DEG_.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  go_files <- list.files(
    src_gsea,
    pattern = "_GO_BP_results\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  re_files <- list.files(
    src_gsea,
    pattern = "_Reactome_results\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  limma_df <- bind_deg_files(deg_files, cluster_col_name)
  go_df    <- bind_gsea_files(go_files, cluster_col_name, "GO:BP")
  re_df    <- bind_gsea_files(re_files, cluster_col_name, "REACTOME")
  
  write_csv_if_present(limma_df, limma_path)
  write_csv_if_present(go_df, go_path)
  write_csv_if_present(re_df, re_path)
}

# =========================
# RUN ALL COMBINATIONS
# =========================
for (analysis_type in analysis_types) {
  for (cluster_col_name in cluster_cols) {
    make_combined_set(analysis_type, cluster_col_name)
  }
}