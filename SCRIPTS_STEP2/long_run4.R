#!/usr/bin/env Rscript
################################################################################
# Comprehensive Single-Cell RNA-seq Analysis Pipeline
# Merges: Setup → QC → Filtering → Annotation → Multi-tool Doublet Detection
# Outputs: Raw object + Final annotated object (2 RDS files only)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(cowplot)
  library(readr)
  library(stringr)
  library(purrr)
  library(patchwork)
  library(scales)
  library(Matrix)
  library(data.table)
  library(SingleR)
  library(celldex, lib="../lib")
  library(BiocParallel)
})

set.seed(1234)

################################################################################
# CONFIGURATION
################################################################################

PROJECT <- "APECED"
BASE_DIR <- "/gpfs/gsfs12/users/cheungf/ACEPD"
OUT_DIR <- file.path(BASE_DIR, "MERGED_PIPELINE_OUTPUT_V5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DONORS <- c("P0001171", "P0001352", "P0001358", "P0002345")
CTRL_RAW <- "P0014164"
CTRL_LAB <- "Control"

QC_THRESHOLDS <- list(
  nCount_min = 1200, nCount_max = 50000,
  nFeature_min = 300, nFeature_max = 7000,
  percent_mito_max = 60,
  percent_ribo_max = 25,
  percent_hb_max = 0.1
)

ANNOTATION_CONFIG <- list(
  p_adj_max = 0.05,
  only_pos = TRUE,
  min_pct = 0.25,
  logfc_threshold = 0.25,
  tolerance = 0.10,
  force_recompute = FALSE
)

################################################################################
# UTILITY FUNCTIONS
################################################################################

log_status <- function(stage, message, level = "INFO") {
  icons <- list("INFO" = "[INFO]", "SUCCESS" = "[OK]", "WARNING" = "[WARN]", 
                "ERROR" = "[ERROR]", "PROGRESS" = "[->]")
  icon <- icons[[level]]
  timestamp <- format(Sys.time(), "%H:%M:%S")
  cat(sprintf("[%s] %s [%s] %s\n", timestamp, icon, stage, message))
}

section_header <- function(title) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("  %s\n", title))
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
}

save_plot_quietly <- function(filename, plot, ...) {
  suppressMessages(ggsave(filename, plot, ...))
  log_status("PLOT", sprintf("Saved: %s", basename(filename)), "SUCCESS")
}

check_and_join_layers <- function(obj, stage_name) {
  current_layers <- SeuratObject::Layers(obj[["RNA"]])
  data_layers <- grep("^data", current_layers, value = TRUE)
  counts_layers <- grep("^counts", current_layers, value = TRUE)
  needs_join <- (length(data_layers) > 1) || (length(counts_layers) > 1)
  if (needs_join) {
    log_status(stage_name, sprintf("Multiple layers detected - joining..."), "INFO")
    if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
      obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
      log_status(stage_name, "Layers joined successfully", "SUCCESS")
    }
  }
  return(obj)
}

join_if_multi <- function(obj, assay = "RNA", stage_name = "JOIN") {
  if (!("JoinLayers" %in% getNamespaceExports("SeuratObject"))) return(obj)
  layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
  data_layers <- grep("^data", layers, value = TRUE)
  counts_layers <- grep("^counts", layers, value = TRUE)
  if (length(data_layers) > 1 || length(counts_layers) > 1) {
    log_status(stage_name, "Collapsing layers...", "PROGRESS")
    obj <- SeuratObject::JoinLayers(obj, assay = assay)
    log_status(stage_name, "Layers joined", "SUCCESS")
  }
  obj
}

################################################################################
# STAGE 1: DATA LOADING & INTEGRATION
################################################################################

stage1_load_and_integrate <- function() {
  section_header("STAGE 1: DATA LOADING & INTEGRATION")
  
  h5_paths <- c(
    "23-301-3_GEX1" = "cellranger_outs/filtered_feature_bc_matrix_23-301-3_GEX1.h5",
    "23-301-3_GEX2" = "cellranger_outs/filtered_feature_bc_matrix_23-301-3_GEX2.h5",
    "23-301-3_GEX3" = "cellranger_outs/filtered_feature_bc_matrix_23-301-3_GEX3.h5",
    "23-301-3_GEX4" = "cellranger_outs/filtered_feature_bc_matrix_23-301-3_GEX4.h5",
    "23-301-3_GEX5" = "cellranger_outs/filtered_feature_bc_matrix_23-301-3_GEX5.h5",
    "23-301-3_GEX6" = "cellranger_outs/filtered_feature_bc_matrix_23-301-3_GEX6.h5"
  )
  
  log_status("LOAD", sprintf("Loading %d GEX samples...", length(h5_paths)), "PROGRESS")
  
  objs <- list()
  for (s in names(h5_paths)) {
    p <- h5_paths[[s]]
    if (!file.exists(p)) {
      log_status("LOAD", sprintf("Missing file: %s", p), "ERROR")
      stop("Missing required file")
    }
    
    log_status("LOAD", sprintf("Reading %s", s), "INFO")
    m <- Read10X_h5(p, use.names = TRUE)
    
    ge <- if (is.list(m) && "Gene Expression" %in% names(m)) {
      m[["Gene Expression"]]
    } else if (inherits(m, "dgCMatrix")) {
      m
    } else {
      stop("Unrecognized Read10X_h5 structure")
    }
    
    clean_bc <- sub("-.*$", "", colnames(ge))
    colnames(ge) <- paste0(clean_bc, "-", s)
    
    sobj <- CreateSeuratObject(counts = ge, project = PROJECT)
    sobj$orig.ident <- s
    sobj <- join_if_multi(sobj, stage_name = paste0("JOIN_", s))
    objs[[s]] <- sobj
    log_status("LOAD", sprintf("%s: %d cells", s, ncol(sobj)), "SUCCESS")
  }
  
  log_status("MERGE", "Merging all samples...", "PROGRESS")
  allsamples <- Reduce(function(a, b) merge(a, b, project = PROJECT), objs)
  allsamples <- join_if_multi(allsamples, stage_name = "JOIN_MERGED")
  
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    allsamples <- SeuratObject::JoinLayers(allsamples)
    log_status("MERGE", "Joined layers (Seurat v5)", "SUCCESS")
  }
  
  log_status("MERGE", sprintf("Total cells: %s", format(ncol(allsamples), big.mark = ",")), "SUCCESS")
  return(allsamples)
}

################################################################################
# STAGE 2: METADATA INTEGRATION (WITH MULTI-TOOL DOUBLETS)
################################################################################

stage2_add_metadata <- function(obj) {
  section_header("STAGE 2: METADATA INTEGRATION")
  
  # ---------------------------
  # QC metrics
  # ---------------------------
  log_status("QC", "Calculating QC metrics...", "PROGRESS")
  obj$percent.mito <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$percent.ribo <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  hb_genes <- grep("^HB", rownames(obj), value = TRUE)
  hb_genes <- hb_genes[!grepl("^HBP", hb_genes)]
  obj$percent.hb <- PercentageFeatureSet(obj, features = hb_genes)
  log_status("QC", "QC metrics calculated", "SUCCESS")
  
  # ---------------------------
  # SNP demux
  # ---------------------------
  log_status("SNP", "Loading SNP demultiplexing data...", "PROGRESS")
  assign_files <- Sys.glob("./DEMUX/*/assignments_refined.tsv")
  
  if (length(assign_files) > 0) {
    read_one <- function(p) {
      smp <- stringr::str_match(p, "(23-301-3_GEX\\d+)_demux")[, 2]
      df <- readr::read_tsv(p, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
      
      # Expect BARCODE + one subject-like column; pick first non-BARCODE
      subj_col <- setdiff(names(df), "BARCODE")[1]
      
      df %>%
        dplyr::transmute(
          core_bc = sub("-.*$", "", BARCODE),
          subject = .data[[subj_col]],
          cell    = paste0(core_bc, "-", smp)
        ) %>%
        dplyr::select(cell, subject)
    }
    
    assign_tbl <- purrr::map_dfr(assign_files, read_one) %>% dplyr::distinct(cell, .keep_all = TRUE)
    
    meta <- obj@meta.data %>% dplyr::mutate(cell = rownames(.))
    
    # Keep order aligned to obj after join
    meta2 <- meta %>%
      dplyr::left_join(assign_tbl, by = "cell", suffix = c("", "_snp")) %>%
      dplyr::mutate(
        Lane = stringr::str_match(orig.ident, "GEX(\\d+)")[, 2],
        SNP_call = dplyr::case_when(
          is.na(subject) ~ "Unassigned",
          stringr::str_detect(as.character(subject), "\\+") ~ "Doublet",
          TRUE ~ "Singleton"
        )
      )
    meta2 <- meta2[match(colnames(obj), meta2$cell), , drop = FALSE]
    
    obj$subject  <- meta2$subject
    obj$Lane     <- meta2$Lane
    obj$SNP_call <- meta2$SNP_call
    log_status("SNP", sprintf("Assigned: %d cells", sum(!is.na(obj$subject))), "SUCCESS")
  } else {
    log_status("SNP", "No SNP assignment files found; skipping.", "WARN")
  }
  
  # ---------------------------
  # Ambient RNA
  # ---------------------------
  log_status("AMBIENT", "Loading ambient RNA scores...", "PROGRESS")
 
  # ---------------------------
  # Ambient RNA  (robust to pre-existing columns)
  # ---------------------------
  log_status("AMBIENT", "Loading ambient RNA scores...", "PROGRESS")
  ambient_files <- Sys.glob("AMBIENT/*GEX*_ambient_scores.tsv")
  
  if (length(ambient_files) > 0) {
    read_ambient <- function(f) {
      smp <- stringr::str_match(basename(f), "(23-301-3_GEX\\d+)_ambient")[, 2]
      df <- readr::read_tsv(f, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
      df %>%
        dplyr::mutate(
          core_bc = sub("-.*$", "", barcode),
          cell = paste0(core_bc, "-", smp),
          ambient_rna_score   = suppressWarnings(as.numeric(ambient_rna_score)),
          percent_top_ambient = suppressWarnings(as.numeric(percent_top_ambient)),
          library_complexity  = suppressWarnings(as.numeric(library_complexity))
        ) %>%
        dplyr::select(cell, ambient_rna_score, percent_top_ambient, library_complexity)
    }
    
    ambient_tbl <- purrr::map_dfr(ambient_files, read_ambient) %>%
      dplyr::distinct(cell, .keep_all = TRUE)
    
    # Use a custom suffix so we never get .x/.y ambiguity
    meta3 <- obj@meta.data %>%
      dplyr::mutate(cell = rownames(.)) %>%
      dplyr::left_join(ambient_tbl, by = "cell", suffix = c("", ".ambient")) %>%
      { .[match(colnames(obj), .$cell), , drop = FALSE] }
    
    # Prefer newly joined values; fall back to any existing ones
    get_col <- function(base) {
      a <- paste0(base, ".ambient")
      if (a %in% names(meta3)) {
        # if an older column exists, coalesce new over old; else just use new
        if (base %in% names(meta3)) {
          dplyr::coalesce(meta3[[a]], meta3[[base]])
        } else {
          meta3[[a]]
        }
      } else if (base %in% names(meta3)) {
        meta3[[base]]
      } else {
        rep(NA_real_, ncol(obj))
      }
    }
    
    obj$ambient_rna_score   <- get_col("ambient_rna_score")
    obj$percent_top_ambient <- get_col("percent_top_ambient")
    obj$library_complexity  <- get_col("library_complexity")
    
    log_status("AMBIENT",
               sprintf("Added scores for %d cells", sum(!is.na(obj$ambient_rna_score))),
               "SUCCESS")
  } else {
    log_status("AMBIENT", "No ambient score files found; skipping.", "WARN")
  }
  
  # ---------------------------
  # MULTI-TOOL DOUBLET DETECTION
  # ---------------------------
  log_status("DOUBLET", "Loading multi-tool doublet detection results...", "PROGRESS")

  
  # Discover input files (prefer per-lane; fallback to ALL)
  doublet_files <- Sys.glob("./DOUBLET/23-301-3_*_combined_DemuxOnly.tsv")
  if (length(doublet_files) == 0) doublet_files <- Sys.glob("./DOUBLET/*combined_DemuxOnly.tsv")
  if (length(doublet_files) == 0) {
    fallback_file <- "/mnt/data/23-301-3_ALL_combined_DemuxOnly.tsv"
    if (file.exists(fallback_file)) doublet_files <- fallback_file
  }
  
  if (length(doublet_files) > 0) {
    
    # ---- Reader that always creates core_bc and normalizes columns
    read_doublets_multi <- function(f) {
      log_status("DOUBLET", sprintf("  Reading %s", basename(f)), "INFO")
      df <- readr::read_tsv(
        f,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        guess_max = 1e6,
        trim_ws = TRUE
      )
      
      # Defensive normalization of the minimal set we need
      if (!"Barcode" %in% names(df)) stop("Doublet TSV missing 'Barcode' column: ", f)
      
      df <- df %>%
        dplyr::mutate(
          Barcode = trimws(Barcode),
          core_bc = toupper(sub("-.*$", "", Barcode)),
          Demuxalot_Individual_Assignment = if ("Demuxalot_Individual_Assignment" %in% names(.)) trimws(Demuxalot_Individual_Assignment) else NA_character_,
          Demuxalot_DropletType           = if ("Demuxalot_DropletType"          %in% names(.)) tolower(trimws(Demuxalot_DropletType))          else NA_character_,
          DoubletDetection_DropletType    = if ("DoubletDetection_DropletType"   %in% names(.)) tolower(trimws(DoubletDetection_DropletType))   else NA_character_,
          scDblFinder_DropletType         = if ("scDblFinder_DropletType"        %in% names(.)) tolower(trimws(scDblFinder_DropletType))        else NA_character_,
          scDblFinder_Score               = if ("scDblFinder_Score"              %in% names(.)) suppressWarnings(as.numeric(scDblFinder_Score)) else NA_real_,
          scds_score                      = if ("scds_score"                     %in% names(.)) suppressWarnings(as.numeric(scds_score))        else NA_real_,
          scds_DropletType                = if ("scds_DropletType"               %in% names(.)) tolower(trimws(scds_DropletType))               else NA_character_,
          scrublet_DropletType            = if ("scrublet_DropletType"           %in% names(.)) tolower(trimws(scrublet_DropletType))           else NA_character_,
          scrublet_Scores                 = if ("scrublet_Scores"                %in% names(.)) suppressWarnings(as.numeric(scrublet_Scores))   else NA_real_
        ) %>%
        dplyr::transmute(
          core_bc,
          Demuxalot_Individual = Demuxalot_Individual_Assignment,
          Demuxalot_Droplet    = Demuxalot_DropletType,
          DoubletDetection_Droplet = DoubletDetection_DropletType,
          scDblFinder_Droplet  = scDblFinder_DropletType,
          scDblFinder_Score,
          scds_Score           = scds_score,
          scds_Droplet         = scds_DropletType,
          scrublet_Droplet     = scrublet_DropletType,
          scrublet_Score       = scrublet_Scores
        )
      
      df
    }
    
    # Read all files and keep a single row per core barcode
    dbl_raw <- purrr::map_dfr(doublet_files, read_doublets_multi) %>%
      dplyr::mutate(core_bc = trimws(core_bc)) %>%
      dplyr::distinct(core_bc, .keep_all = TRUE)
    
    # ---- Build lane candidates from the Seurat object itself
    obj_names_uc <- toupper(trimws(colnames(obj)))
    lanes_in_obj <- unique(sub("^[^-]+-", "", obj_names_uc))  # e.g., "23-301-3_GEX1"
    
    # For each lane, see if "<core>-<lane>" exists in the object
    cand_mat <- if (length(lanes_in_obj) > 0 && nrow(dbl_raw) > 0) {
      sapply(lanes_in_obj, function(ln) paste0(dbl_raw$core_bc, "-", ln) %in% obj_names_uc)
    } else {
      matrix(FALSE, nrow = nrow(dbl_raw), ncol = 0)
    }
    
    # Rows that hit at least one lane; choose first hit per row
    any_hit <- if (ncol(cand_mat) > 0) rowSums(cand_mat) > 0 else rep(FALSE, nrow(dbl_raw))
    chosen_idx  <- if (ncol(cand_mat) > 0) max.col(cand_mat, ties.method = "first") else integer(nrow(dbl_raw))
    chosen_lane <- if (length(lanes_in_obj) > 0) lanes_in_obj[chosen_idx] else rep(NA_character_, nrow(dbl_raw))
    chosen_lane[!any_hit] <- NA_character_
    
    placed <- any_hit
    n_matched <- sum(placed)
    
    if (n_matched == 0) {
      log_status("DOUBLET", "Matched 0 cells after trimming; barcodes may not overlap this object.", "WARN")
    }
    
    # Keep only rows we can place and build exact Seurat key
    dbl_tbl <- dbl_raw[placed, , drop = FALSE] %>%
      dplyr::mutate(cell = paste0(core_bc, "-", chosen_lane[placed])) %>%
      dplyr::select(-core_bc)
    
    # ---- Consensus across tools
    droplet_cols <- intersect(
      c("Demuxalot_Droplet","DoubletDetection_Droplet","scDblFinder_Droplet","scds_Droplet","scrublet_Droplet"),
      colnames(dbl_tbl)
    )
    
    if (length(droplet_cols) > 0) {
      dbl_mat    <- as.matrix(dbl_tbl[droplet_cols] == "doublet")
      called_mat <- as.matrix(!is.na(dbl_tbl[droplet_cols]))
      n_doublet_calls <- rowSums(dbl_mat, na.rm = TRUE)
      n_tools_called  <- rowSums(called_mat, na.rm = TRUE)
    } else {
      n_doublet_calls <- integer(nrow(dbl_tbl))
      n_tools_called  <- integer(nrow(dbl_tbl))
    }
    
    dbl_tbl <- dbl_tbl %>%
      dplyr::mutate(
        n_doublet_calls = n_doublet_calls,
        n_tools_called  = n_tools_called,
        doublet_consensus = dplyr::case_when(
          n_tools_called == 0 ~ "unknown",
          n_doublet_calls / pmax(n_tools_called, 1L) >= 0.5 ~ "doublet",
          TRUE ~ "singlet"
        ),
        doublet_confidence = dplyr::if_else(
          n_tools_called > 0,
          pmax(n_doublet_calls, n_tools_called - n_doublet_calls) / n_tools_called,
          as.numeric(NA)
        )
      )
    
    # ---- Join back with predictable suffix and coalesce
    meta_doublet <- obj@meta.data %>%
      dplyr::mutate(cell = rownames(.)) %>%
      dplyr::left_join(dbl_tbl, by = "cell", suffix = c("", ".dbl")) %>%
      { .[match(colnames(obj), .$cell), , drop = FALSE] }
    
    get_dbl <- function(base) {
      new <- paste0(base, ".dbl")
      if (new %in% names(meta_doublet) || base %in% names(meta_doublet)) {
        dplyr::coalesce(meta_doublet[[new]], meta_doublet[[base]])
      } else {
        rep(NA, ncol(obj))
      }
    }
    
    obj$Demuxalot_Individual     <- get_dbl("Demuxalot_Individual")
    obj$Demuxalot_Droplet        <- get_dbl("Demuxalot_Droplet")
    obj$DoubletDetection_Droplet <- get_dbl("DoubletDetection_Droplet")
    obj$scDblFinder_Droplet      <- get_dbl("scDblFinder_Droplet")
    obj$scDblFinder_Score        <- get_dbl("scDblFinder_Score")
    obj$scds_Score               <- get_dbl("scds_Score")
    obj$scds_Droplet             <- get_dbl("scds_Droplet")
    obj$scrublet_Droplet         <- get_dbl("scrublet_Droplet")
    obj$scrublet_Score           <- get_dbl("scrublet_Score")
    obj$doublet_consensus        <- get_dbl("doublet_consensus")
    obj$doublet_confidence       <- get_dbl("doublet_confidence")
    obj$n_doublet_calls          <- get_dbl("n_doublet_calls")
    obj$n_tools_called           <- get_dbl("n_tools_called")
    
    log_status(
      "DOUBLET",
      sprintf("Attached doublet calls to %s / %s cells (cores placed: %s)",
              format(sum(!is.na(obj$Demuxalot_Droplet) | !is.na(obj$DoubletDetection_Droplet) |
                           !is.na(obj$scDblFinder_Droplet) | !is.na(obj$scds_Droplet) |
                           !is.na(obj$scrublet_Droplet)), big.mark = ","),
              format(ncol(obj), big.mark = ","),
              format(n_matched, big.mark = ",")
      ),
      if (n_matched > 0) "SUCCESS" else "WARN"
    )
    
    cat("\n")
    for (tool_col in c("DoubletDetection_Droplet", "scDblFinder_Droplet", "scds_Droplet", "scrublet_Droplet")) {
      if (tool_col %in% colnames(obj@meta.data)) {
        tool_name <- sub("_Droplet$", "", tool_col)
        denom <- sum(!is.na(obj@meta.data[[tool_col]]))
        if (denom > 0) {
          n_dbl <- sum(obj@meta.data[[tool_col]] == "doublet", na.rm = TRUE)
          pct_dbl <- 100 * n_dbl / denom
          log_status("DOUBLET", sprintf("  %s: %.1f%%", tool_name, pct_dbl), "INFO")
        } else {
          log_status("DOUBLET", sprintf("  %s: no calls", tool_name), "INFO")
        }
      }
    }
    cat("\n")
    
  } else {
    log_status("DOUBLET", "No doublet files found; skipping.", "WARN")
  }
  return(obj)
  }

################################################################################
# STAGE 3: SINGLETON FILTERING & HARMONIZATION
################################################################################

stage3_filter_and_harmonize <- function(obj) {
  section_header("STAGE 3: SINGLETON FILTERING & HARMONIZATION")
  
  md <- obj@meta.data
  md$Subject_raw <- as.character(md$subject)
  md$Subject_fixed <- md$Subject_raw
  md$Subject_fixed[md$Subject_fixed == CTRL_RAW] <- CTRL_LAB
  
  is_singleton <- tolower(as.character(md$SNP_call)) == "singleton"
  no_plus <- !grepl("\\+", md$Subject_raw, perl = TRUE)
  in_whitelist <- md$Subject_fixed %in% c(DONORS, CTRL_LAB)
  keep <- is_singleton & no_plus & in_whitelist
  
  log_status("FILTER", sprintf("Singleton filter: %d/%d cells (%.1f%%)", 
                               sum(keep), nrow(md), 100 * mean(keep)), "INFO")
  
  obj@meta.data <- md
  obj <- subset(obj, cells = rownames(obj@meta.data)[keep])
  obj <- join_if_multi(obj, stage_name = "JOIN_AFTER_SINGLETON_SUBSET")
  md <- obj@meta.data
  
  lane_tp_map <- function(lane, subj) {
    lane <- as.character(lane)
    subj <- as.character(subj)
    if (is.na(lane) || is.na(subj)) return("Unknown")
    if (subj == CTRL_LAB) return("Spike")
    if (!(subj %in% DONORS)) return("Unknown")
    odd <- lane %in% c("1", "3", "5")
    if (odd) {
      if (subj %in% c("P0001171", "P0001352")) "Before" else "After"
    } else {
      if (subj %in% c("P0001171", "P0001352")) "After" else "Before"
    }
  }
  
  md$Lane <- as.character(md$Lane)
  md$TimePoint_fixed <- mapply(lane_tp_map, md$Lane, md$Subject_fixed, USE.NAMES = FALSE)
  md$Treatment <- case_when(
    md$TimePoint_fixed == "Before" ~ "Untreated",
    md$TimePoint_fixed == "After" ~ "Treated",
    md$TimePoint_fixed == "Spike" ~ "Spike",
    TRUE ~ "Unknown")
  
  obj@meta.data <- md
  log_status("FILTER", sprintf("Final: %s cells", format(ncol(obj), big.mark = ",")), "SUCCESS")
  return(obj)
}

################################################################################
# STAGE 4: DIMENSIONAL REDUCTION
################################################################################

stage4_dimension_reduction <- function(obj) {
  section_header("STAGE 4: DIMENSIONAL REDUCTION")
  
  obj <- join_if_multi(obj, stage_name = "JOIN_BEFORE_DIMRED")
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  }
  
  log_status("DIMRED", "Normalizing...", "PROGRESS")
  obj <- NormalizeData(obj, verbose = FALSE)
  
  log_status("DIMRED", "Finding variable features...", "PROGRESS")
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  
  log_status("DIMRED", "Scaling...", "PROGRESS")
  obj <- ScaleData(obj, verbose = FALSE)
  
  log_status("DIMRED", "Running PCA...", "PROGRESS")
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  
  log_status("DIMRED", "Computing neighbors and clusters...", "PROGRESS")
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  log_status("DIMRED", "Running UMAP...", "PROGRESS")
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  
  log_status("DIMRED", sprintf("Identified %d clusters", length(unique(obj$seurat_clusters))), "SUCCESS")
  return(obj)
}

################################################################################
# STAGE 5: QC FILTERING
################################################################################

stage5_qc_filtering <- function(obj) {
  section_header("STAGE 5: QC FILTERING")
  
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    obj$percent.mt <- obj$percent.mito
  }
  
  if (!"percent.ribo" %in% colnames(obj@meta.data)) {
    ribo <- grep("^RPL|^RPS", rownames(obj), perl = TRUE, value = TRUE)
    if (length(ribo) > 0) {
      counts_mat <- GetAssayData(obj, slot = "counts")
      ribo_counts <- Matrix::colSums(counts_mat[ribo, , drop = FALSE])
      total_counts <- Matrix::colSums(counts_mat)
      obj$percent.ribo <- (ribo_counts / pmax(total_counts, 1)) * 100
    }
  }
  
  pre_filter <- ncol(obj)
  obj_filtered <- subset(obj,
                         subset = nCount_RNA > QC_THRESHOLDS$nCount_min & nCount_RNA < QC_THRESHOLDS$nCount_max &
                           nFeature_RNA > QC_THRESHOLDS$nFeature_min & nFeature_RNA < QC_THRESHOLDS$nFeature_max &
                           percent.mt < QC_THRESHOLDS$percent_mito_max &
                           percent.ribo < QC_THRESHOLDS$percent_ribo_max &
                           percent.hb < QC_THRESHOLDS$percent_hb_max)
  
  obj_filtered <- join_if_multi(obj_filtered, stage_name = "JOIN_AFTER_QC")
  post_filter <- ncol(obj_filtered)
  
  log_status("FILTER", sprintf("%s/%s cells retained (%.1f%%)", 
                               format(post_filter, big.mark = ","),
                               format(pre_filter, big.mark = ","),
                               100 * post_filter / pre_filter), "SUCCESS")
  return(obj_filtered)
}

################################################################################
# STAGE 6: CD45 GATING
################################################################################

stage6_cd45_gating <- function(obj) {
  section_header("STAGE 6: CD45/EPCAM GATING")
  
  rn <- rownames(obj)
  cd45_gene <- if ("PTPRC" %in% rn) "PTPRC" else if ("CD45" %in% rn) "CD45" else NA_character_
  epcam_gene <- if ("EPCAM" %in% rn) "EPCAM" else NA_character_
  
  obj$CD45_expr <- if (!is.na(cd45_gene)) FetchData(obj, cd45_gene)[, 1] else 0
  obj$EPCAM_expr <- if (!is.na(epcam_gene)) FetchData(obj, epcam_gene)[, 1] else 0
  
  thr_kmeans <- function(x) {
    x0 <- x[x > 0]
    if (length(unique(x0)) < 2) return(0)
    km <- kmeans(x0, centers = 2, nstart = 10)
    m <- tapply(x0, km$cluster, mean)
    mean(sort(m))
  }
  
  thr_cd45 <- thr_kmeans(obj$CD45_expr)
  thr_epcam <- thr_kmeans(obj$EPCAM_expr)
  obj$CD45pos <- (obj$CD45_expr > thr_cd45) & (obj$EPCAM_expr <= thr_epcam)
  
  log_status("GATE", sprintf("CD45+ cells: %s (%.1f%%)", 
                             format(sum(obj$CD45pos), big.mark = ","), 
                             100 * sum(obj$CD45pos) / ncol(obj)), "SUCCESS")
  return(obj)
}

################################################################################
# STAGE 7: SINGLER ANNOTATION
################################################################################

stage7_singler_annotation <- function(obj) {
  section_header("STAGE 7: SINGLER REFERENCE ANNOTATION")
  
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    log_status("SINGLER", "SingleR not installed - skipping reference annotation", "WARNING")
    return(obj)
  }
  
  if (!requireNamespace("celldex", quietly = TRUE)) {
    log_status("SINGLER", "celldex not installed - skipping reference annotation", "WARNING")
    return(obj)
  }
  
  obj <- check_and_join_layers(obj, "SINGLER")
  
  log_status("SINGLER", "Extracting normalized expression matrix...", "PROGRESS")
  
  data_layers <- grep("^data", SeuratObject::Layers(obj[["RNA"]]), value = TRUE)
  
  if (length(data_layers) > 0) {
    counts_for_singler <- SeuratObject::LayerData(obj[["RNA"]], layer = data_layers[1])
  } else {
    log_status("SINGLER", "No normalized data found - normalizing now", "INFO")
    obj <- NormalizeData(obj, verbose = FALSE)
    data_layers <- grep("^data", SeuratObject::Layers(obj[["RNA"]]), value = TRUE)
    counts_for_singler <- SeuratObject::LayerData(obj[["RNA"]], layer = data_layers[1])
  }
  
  log_status("SINGLER", sprintf("Matrix: %d genes × %d cells", 
                                nrow(counts_for_singler), 
                                ncol(counts_for_singler)), "INFO")
  
  # Monaco Reference
  log_status("SINGLER", "Loading Monaco Immune reference...", "PROGRESS")
  
  tryCatch({
    ref_monaco <- celldex::MonacoImmuneData()
    log_status("SINGLER", sprintf("Monaco reference loaded: %d cell types", 
                                  length(unique(ref_monaco$label.main))), "SUCCESS")
    
    log_status("SINGLER", "Running SingleR with Monaco (main labels)...", "PROGRESS")
    pred_monaco_main <- SingleR(
      test = counts_for_singler,
      ref = ref_monaco,
      labels = ref_monaco$label.main,
      BPPARAM = BiocParallel::SerialParam()
    )
    
    log_status("SINGLER", "Running SingleR with Monaco (fine labels)...", "PROGRESS")
    pred_monaco_fine <- SingleR(
      test = counts_for_singler,
      ref = ref_monaco,
      labels = ref_monaco$label.fine,
      BPPARAM = BiocParallel::SerialParam()
    )
    
    obj$Monaco_main <- pred_monaco_main$labels
    obj$Monaco_main_pruned <- pred_monaco_main$pruned.labels
    obj$Monaco_main_scores <- pred_monaco_main$scores[cbind(seq_len(nrow(pred_monaco_main)), 
                                                            match(pred_monaco_main$labels, colnames(pred_monaco_main$scores)))]
    
    obj$Monaco_fine <- pred_monaco_fine$labels
    obj$Monaco_fine_pruned <- pred_monaco_fine$pruned.labels
    obj$Monaco_fine_scores <- pred_monaco_fine$scores[cbind(seq_len(nrow(pred_monaco_fine)), 
                                                            match(pred_monaco_fine$labels, colnames(pred_monaco_fine$scores)))]
    
    n_main <- length(unique(obj$Monaco_main))
    n_pruned_main <- sum(!is.na(obj$Monaco_main_pruned))
    pct_pruned_main <- 100 * n_pruned_main / ncol(obj)
    
    log_status("SINGLER", sprintf("Monaco main: %d unique labels, %.1f%% cells retained after pruning",
                                  n_main, pct_pruned_main), "SUCCESS")
    
    n_fine <- length(unique(obj$Monaco_fine))
    n_pruned_fine <- sum(!is.na(obj$Monaco_fine_pruned))
    pct_pruned_fine <- 100 * n_pruned_fine / ncol(obj)
    
    log_status("SINGLER", sprintf("Monaco fine: %d unique labels, %.1f%% cells retained after pruning",
                                  n_fine, pct_pruned_fine), "SUCCESS")
    
  }, error = function(e) {
    log_status("SINGLER", sprintf("Monaco annotation failed: %s", e$message), "WARNING")
  })
  
  # HPCA Reference
  log_status("SINGLER", "Loading Human Primary Cell Atlas (HPCA) reference...", "PROGRESS")
  
  tryCatch({
    ref_hpca <- celldex::HumanPrimaryCellAtlasData()
    log_status("SINGLER", sprintf("HPCA reference loaded: %d cell types", 
                                  length(unique(ref_hpca$label.main))), "SUCCESS")
    
    log_status("SINGLER", "Running SingleR with HPCA (main labels)...", "PROGRESS")
    pred_hpca_main <- SingleR(
      test = counts_for_singler,
      ref = ref_hpca,
      labels = ref_hpca$label.main,
      BPPARAM = BiocParallel::SerialParam()
    )
    
    log_status("SINGLER", "Running SingleR with HPCA (fine labels)...", "PROGRESS")
    pred_hpca_fine <- SingleR(
      test = counts_for_singler,
      ref = ref_hpca,
      labels = ref_hpca$label.fine,
      BPPARAM = BiocParallel::SerialParam()
    )
    
    obj$HPCA_main <- pred_hpca_main$labels
    obj$HPCA_main_pruned <- pred_hpca_main$pruned.labels
    obj$HPCA_main_scores <- pred_hpca_main$scores[cbind(seq_len(nrow(pred_hpca_main)), 
                                                        match(pred_hpca_main$labels, colnames(pred_hpca_main$scores)))]
    
    obj$HPCA_fine <- pred_hpca_fine$labels
    obj$HPCA_fine_pruned <- pred_hpca_fine$pruned.labels
    obj$HPCA_fine_scores <- pred_hpca_fine$scores[cbind(seq_len(nrow(pred_hpca_fine)), 
                                                        match(pred_hpca_fine$labels, colnames(pred_hpca_fine$scores)))]
    
    n_main <- length(unique(obj$HPCA_main))
    n_pruned_main <- sum(!is.na(obj$HPCA_main_pruned))
    pct_pruned_main <- 100 * n_pruned_main / ncol(obj)
    
    log_status("SINGLER", sprintf("HPCA main: %d unique labels, %.1f%% cells retained after pruning",
                                  n_main, pct_pruned_main), "SUCCESS")
    
    n_fine <- length(unique(obj$HPCA_fine))
    n_pruned_fine <- sum(!is.na(obj$HPCA_fine_pruned))
    pct_pruned_fine <- 100 * n_pruned_fine / ncol(obj)
    
    log_status("SINGLER", sprintf("HPCA fine: %d unique labels, %.1f%% cells retained after pruning",
                                  n_fine, pct_pruned_fine), "SUCCESS")
    
  }, error = function(e) {
    log_status("SINGLER", sprintf("HPCA annotation failed: %s", e$message), "WARNING")
  })
  
  # Save SingleR results
  singler_dir <- file.path(OUT_DIR, "singler_results")
  dir.create(singler_dir, showWarnings = FALSE)
  
  singler_cols <- grep("^(Monaco|HPCA)_", colnames(obj@meta.data), value = TRUE)
  if (length(singler_cols) > 0) {
    singler_df <- obj@meta.data[, singler_cols, drop = FALSE]
    singler_df$cell_id <- rownames(singler_df)
    fwrite(singler_df, file.path(singler_dir, "singler_per_cell_annotations.csv"))
    log_status("SINGLER", "SAVED: singler_per_cell_annotations.csv", "SUCCESS")
  }
  
  return(obj)
}

################################################################################
# STAGE 8: CLUSTER ANNOTATION WITH SMART CACHING
################################################################################

stage8_cluster_annotation <- function(obj) {
  section_header("STAGE 8: CLUSTER ANNOTATION")
  
  markers_dir <- file.path(OUT_DIR, "markers")
  dir.create(markers_dir, showWarnings = FALSE, recursive = TRUE)
  
  markers_all_file <- file.path(markers_dir, "FindAllMarkers_all_results.csv")
  markers_top5_file <- file.path(markers_dir, "FindAllMarkers_top5_per_cluster.csv")
  markers_sig_file <- file.path(markers_dir, "FindAllMarkers_significant_only.csv")
  
  obj <- join_if_multi(obj, stage_name = "JOIN_BEFORE_CLUSTER_ANNOT")
  obj <- check_and_join_layers(obj, "ANNOT")
  
  cluster_col <- "seurat_clusters"
  obj$cluster <- factor(obj@meta.data[[cluster_col]])
  
  log_status("ANNOT", "Computing cluster statistics...", "PROGRESS")
  
  epi_genes <- intersect(c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"), rownames(obj))
  if (length(epi_genes) > 0 && !"EpithelialScore" %in% colnames(obj@meta.data)) {
    obj <- AddModuleScore(obj, features = list(epi_genes), name = "EpithelialScore")
    obj$EpithelialScore <- obj$EpithelialScore1
    obj$EpithelialScore1 <- NULL
  } else if (!"EpithelialScore" %in% colnames(obj@meta.data)) {
    obj$EpithelialScore <- 0
  }
  
  if ("percent.mt" %in% colnames(obj@meta.data)) {
    pct_mt <- obj@meta.data$percent.mt
  } else if ("percent.mito" %in% colnames(obj@meta.data)) {
    pct_mt <- obj@meta.data$percent.mito
  } else {
    pct_mt <- rep(0, ncol(obj))
  }
  
  md <- obj@meta.data
  md$cluster <- as.character(obj$cluster)
  md$percent_mt_safe <- pct_mt
  
  stats_df <- md %>%
    group_by(cluster) %>%
    summarise(
      n_cells = n(),
      median_mt_percent = median(percent_mt_safe, na.rm = TRUE),
      mean_mt_percent = mean(percent_mt_safe, na.rm = TRUE),
      median_epi_score = median(EpithelialScore, na.rm = TRUE),
      mean_epi_score = mean(EpithelialScore, na.rm = TRUE),
      median_cd45_expr = median(CD45_expr, na.rm = TRUE),
      mean_cd45_expr = mean(CD45_expr, na.rm = TRUE),
      .groups = "drop"
    )
  
  decide_lineage <- function(med_cd45, med_epi, tol = ANNOTATION_CONFIG$tolerance) {
    if (is.na(med_cd45) || is.na(med_epi)) return(c("Mixed", "NA"))
    if (med_cd45 == 0) return(c("Non-lymphoid", "med_cd45==0"))
    if (med_cd45 > med_epi + tol) return(c("Lymphoid", "med_cd45>med_epi+tol"))
    if (med_cd45 < med_epi - tol) return(c("Non-lymphoid", "med_cd45<med_epi-tol"))
    c("Mixed", "within_tolerance")
  }
  
  lineage_res <- mapply(decide_lineage, stats_df$median_cd45_expr, stats_df$median_epi_score)
  stats_df$lineage_call <- lineage_res[1, ]
  stats_df$lineage_rule <- lineage_res[2, ]
  stats_df$pct_of_total <- round(100 * stats_df$n_cells / ncol(obj), 2)
  stats_df$high_mt_flag <- ifelse(stats_df$median_mt_percent > 20, "HIGH", "OK")
  
  log_status("ANNOT", sprintf("Lineage: %d Lymphoid, %d Non-lymphoid, %d Mixed",
                              sum(stats_df$lineage_call == "Lymphoid"),
                              sum(stats_df$lineage_call == "Non-lymphoid"),
                              sum(stats_df$lineage_call == "Mixed")), "SUCCESS")
  
  log_status("ANNOT", "Checking for existing marker results...", "PROGRESS")
  
  load_cached <- !ANNOTATION_CONFIG$force_recompute && 
    file.exists(markers_all_file) &&
    file.exists(markers_top5_file)
  
  if (load_cached) {
    log_status("ANNOT", "[LOAD] Loading cached FindAllMarkers results...", "INFO")
    all_markers <- fread(markers_all_file)
    markers_df <- fread(markers_top5_file)
    log_status("ANNOT", sprintf("[OK] Loaded %s markers from cache", format(nrow(all_markers), big.mark=",")), "SUCCESS")
  } else {
    log_status("ANNOT", "Finding marker genes (this may take several minutes)...", "PROGRESS")
    
    all_markers <- tryCatch({
      FindAllMarkers(
        object = obj,
        only.pos = ANNOTATION_CONFIG$only_pos,
        min.pct = ANNOTATION_CONFIG$min_pct,
        logfc.threshold = ANNOTATION_CONFIG$logfc_threshold,
        return.thresh = ANNOTATION_CONFIG$p_adj_max,
        verbose = FALSE
      )
    }, error = function(e) {
      log_status("ANNOT", sprintf("Marker finding error: %s", e$message), "WARNING")
      NULL
    })
    
    if (!is.null(all_markers) && nrow(all_markers) > 0) {
      all_markers$cluster <- as.character(all_markers$cluster)
      
      fwrite(all_markers, markers_all_file)
      log_status("ANNOT", sprintf("SAVED: %s (%s markers)", 
                                  basename(markers_all_file), 
                                  format(nrow(all_markers), big.mark=",")), "SUCCESS")
      
      sig_markers <- all_markers %>% filter(p_val_adj <= ANNOTATION_CONFIG$p_adj_max)
      fwrite(sig_markers, markers_sig_file)
      log_status("ANNOT", sprintf("SAVED: %s (%s significant)", 
                                  basename(markers_sig_file), 
                                  format(nrow(sig_markers), big.mark=",")), "SUCCESS")
      
      markers_df <- all_markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        slice_head(n = 5) %>%
        summarise(
          top_markers = paste0(gene, collapse = ", "),
          top_logFC = paste0(round(avg_log2FC, 2), collapse = ", "),
          .groups = "drop"
        )
      
      fwrite(markers_df, markers_top5_file)
      log_status("ANNOT", sprintf("SAVED: %s", basename(markers_top5_file)), "SUCCESS")
      
      log_status("ANNOT", "Saving per-cluster marker files...", "PROGRESS")
      for (clust in unique(all_markers$cluster)) {
        clust_markers <- all_markers %>% filter(cluster == clust)
        clust_file <- file.path(markers_dir, sprintf("cluster_%s_markers.csv", clust))
        fwrite(clust_markers, clust_file)
      }
      log_status("ANNOT", sprintf("SAVED: Per-cluster files for %d clusters", 
                                  length(unique(all_markers$cluster))), "SUCCESS")
      
    } else {
      markers_df <- data.frame(cluster = character(), top_markers = character(), 
                               top_logFC = character())
      log_status("ANNOT", "No significant markers found", "WARNING")
    }
  }
  
  log_status("ANNOT", "Extracting reference annotations...", "PROGRESS")
  
  monaco_cols <- c("Monaco_main", "Monaco_main_pruned", "Monaco_fine", "Monaco_fine_pruned")
  available_monaco <- monaco_cols[monaco_cols %in% colnames(md)]
  
  if (length(available_monaco) > 0) {
    get_mode <- function(x) {
      ux <- unique(x[!is.na(x)])
      if (length(ux) == 0) return(NA_character_)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    monaco_df <- md %>%
      group_by(cluster) %>%
      summarise(
        monaco_main = if("Monaco_main" %in% available_monaco) get_mode(Monaco_main) else NA_character_,
        monaco_main_pruned = if("Monaco_main_pruned" %in% available_monaco) get_mode(Monaco_main_pruned) else NA_character_,
        .groups = "drop"
      )
    log_status("ANNOT", "Extracted Monaco annotations", "SUCCESS")
  } else {
    monaco_df <- data.frame(cluster = unique(as.character(md$cluster)),
                            monaco_main = NA_character_,
                            monaco_main_pruned = NA_character_)
  }
  
  hpca_cols <- c("HPCA_main", "HPCA_main_pruned", "HPCA_fine", "HPCA_fine_pruned")
  available_hpca <- hpca_cols[hpca_cols %in% colnames(md)]
  
  if (length(available_hpca) > 0) {
    hpca_df <- md %>%
      group_by(cluster) %>%
      summarise(
        hpca_main = if("HPCA_main" %in% available_hpca) get_mode(HPCA_main) else NA_character_,
        hpca_main_pruned = if("HPCA_main_pruned" %in% available_hpca) get_mode(HPCA_main_pruned) else NA_character_,
        .groups = "drop"
      )
    log_status("ANNOT", "Extracted HPCA annotations", "SUCCESS")
  } else {
    hpca_df <- data.frame(cluster = unique(as.character(md$cluster)),
                          hpca_main = NA_character_,
                          hpca_main_pruned = NA_character_)
  }
  
  stats_df$cluster <- as.character(stats_df$cluster)
  markers_df$cluster <- as.character(markers_df$cluster)
  monaco_df$cluster <- as.character(monaco_df$cluster)
  hpca_df$cluster <- as.character(hpca_df$cluster)
  
  final_annot <- stats_df %>%
    left_join(markers_df, by = "cluster") %>%
    left_join(monaco_df, by = "cluster") %>%
    left_join(hpca_df, by = "cluster")
  
  final_annot$annotation_summary <- with(final_annot, paste0(
    "Lineage: ", lineage_call,
    " | Monaco: ", ifelse(is.na(monaco_main), "NA", monaco_main),
    " | HPCA: ", ifelse(is.na(hpca_main), "NA", hpca_main)
  ))
  
  annot_file <- file.path(OUT_DIR, "cluster_annotations_comprehensive.csv")
  fwrite(final_annot, annot_file)
  log_status("ANNOT", sprintf("SAVED: %s", basename(annot_file)), "SUCCESS")
  
  annot_lookup <- setNames(final_annot$annotation_summary, final_annot$cluster)
  current_clusters <- as.character(obj@meta.data$cluster)
  obj@meta.data$cluster_annotation <- annot_lookup[current_clusters]
  
  if (!"cluster_annotation" %in% colnames(obj@meta.data)) {
    log_status("ANNOT", "Failed to add cluster_annotation to metadata", "ERROR")
  } else {
    n_annotated <- sum(!is.na(obj@meta.data$cluster_annotation))
    log_status("ANNOT", sprintf("[OK] Added cluster_annotation: %d/%d cells", 
                                n_annotated, ncol(obj)), "SUCCESS")
  }
  
  return(obj)
}

################################################################################
# STAGE 9: VISUALIZATION (UPDATED WITH DOUBLET PLOTS)
################################################################################

stage9_generate_plots <- function(obj) {
  section_header("STAGE 9: GENERATING VISUALIZATIONS")
  
  plots_dir <- file.path(OUT_DIR, "plots")
  dir.create(plots_dir, showWarnings = FALSE)
  
  doublet_plots_dir <- file.path(plots_dir, "doublet_detection")
  dir.create(doublet_plots_dir, showWarnings = FALSE)
  
  subject_pal <- c(
    "P0001171" = "#F06565",
    "P0001352" = "#9ACD32",
    "P0001358" = "#00B5D8",
    "P0002345" = "#B085F5",
    "Control" = "#333333"
  )
  
  log_status("PLOT", "Generating UMAP plots...", "PROGRESS")
  
  p1 <- DimPlot(obj, reduction = "umap", group.by = "Subject_fixed",
                cols = subject_pal, pt.size = 0.2, raster = TRUE) +
    ggtitle("UMAP colored by Subject")
  save_plot_quietly(file.path(plots_dir, "UMAP_by_Subject.png"), p1, 
                    width = 8, height = 6.5, dpi = 300)
  
  if ("TimePoint_fixed" %in% colnames(obj@meta.data)) {
    p2 <- DimPlot(obj, reduction = "umap", group.by = "TimePoint_fixed",
                  pt.size = 0.2, raster = TRUE) +
      ggtitle("UMAP colored by TimePoint")
    save_plot_quietly(file.path(plots_dir, "UMAP_by_TimePoint.png"), p2,
                      width = 8, height = 6.5, dpi = 300)
  }
  
  p3 <- DimPlot(obj, reduction = "umap", group.by = "cluster",
                pt.size = 0.2, raster = TRUE, label = TRUE) +
    ggtitle("UMAP colored by Cluster")
  save_plot_quietly(file.path(plots_dir, "UMAP_by_Cluster.png"), p3,
                    width = 9, height = 7, dpi = 300)
  
  p4 <- DimPlot(obj, reduction = "umap", group.by = "Subject_fixed",
                split.by = "Lane", cols = subject_pal, ncol = 3,
                pt.size = 0.2, raster = TRUE) +
    ggtitle("UMAP by Subject, split by Lane")
  save_plot_quietly(file.path(plots_dir, "UMAP_by_Subject_split_Lane.png"), p4,
                    width = 16, height = 10, dpi = 300)
  
  if ("CD45pos" %in% colnames(obj@meta.data)) {
    p5 <- DimPlot(obj, reduction = "umap", group.by = "CD45pos",
                  pt.size = 0.2, raster = TRUE) +
      ggtitle("UMAP: CD45+ gate (PTPRC high, EPCAM low)")
    save_plot_quietly(file.path(plots_dir, "UMAP_CD45_gate.png"), p5,
                      width = 8, height = 6.5, dpi = 300)
  }
  
  # ========== DOUBLET DETECTION PLOTS ==========
  log_status("PLOT", "Generating doublet detection plots...", "PROGRESS")
  
  doublet_colors <- c("singlet" = "#2E7D32", "doublet" = "#D32F2F", "unknown" = "#757575")
  
  if ("doublet_consensus" %in% colnames(obj@meta.data)) {
    p_consensus <- DimPlot(obj, reduction = "umap", group.by = "doublet_consensus",
                           cols = doublet_colors, pt.size = 1, raster = TRUE) +
      ggtitle("UMAP: Consensus Doublet Detection") +
      theme(legend.position = "right")
    save_plot_quietly(file.path(doublet_plots_dir, "UMAP_doublet_consensus.png"), 
                      p_consensus, width = 9, height = 7, dpi = 300)
  }
  
  doublet_tools <- c("DoubletDetection_Droplet", "scDblFinder_Droplet", 
                     "scds_Droplet", "scrublet_Droplet")
  tool_names <- c("DoubletDetection", "scDblFinder", "scds", "scrublet")
  
  available_tools <- doublet_tools[doublet_tools %in% colnames(obj@meta.data)]
  
  if (length(available_tools) > 0) {
    plot_list <- list()
    
    for (i in seq_along(available_tools)) {
      tool_col <- available_tools[i]
      tool_name <- tool_names[i]
      
      p_tool <- DimPlot(obj, reduction = "umap", group.by = tool_col,
                        cols = doublet_colors, pt.size = 1, raster = TRUE) +
        ggtitle(tool_name) +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.5, face = "bold"))
      
      plot_list[[tool_name]] <- p_tool
      
      save_plot_quietly(
        file.path(doublet_plots_dir, sprintf("UMAP_doublet_%s.png", tool_name)),
        p_tool, width = 7, height = 6, dpi = 300
      )
    }
    
    if (length(plot_list) == 4) {
      p_combined <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(
          title = "Doublet Detection: Comparison Across Tools",
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        )
      save_plot_quietly(file.path(doublet_plots_dir, "UMAP_doublet_all_tools_comparison.png"),
                        p_combined, width = 14, height = 12, dpi = 300)
    }
  }
  
  if ("doublet_confidence" %in% colnames(obj@meta.data)) {
    p_confidence <- FeaturePlot(obj, features = "doublet_confidence", 
                                reduction = "umap", pt.size = 1, raster = TRUE) +
      scale_color_gradient2(low = "#D32F2F", mid = "#FDD835", high = "#2E7D32",
                            midpoint = 0.5, name = "Confidence") +
      ggtitle("UMAP: Doublet Detection Confidence Score") +
      theme(legend.position = "right")
    save_plot_quietly(file.path(doublet_plots_dir, "UMAP_doublet_confidence.png"),
                      p_confidence, width = 9, height = 7, dpi = 300)
  }
  
  if ("n_doublet_calls" %in% colnames(obj@meta.data)) {
    p_n_calls <- FeaturePlot(obj, features = "n_doublet_calls",
                             reduction = "umap", pt.size = 1, raster = TRUE) +
      scale_color_viridis_c(option = "plasma", name = "# Tools\nCalling\nDoublet") +
      ggtitle("UMAP: Number of Tools Calling Doublet") +
      theme(legend.position = "right")
    save_plot_quietly(file.path(doublet_plots_dir, "UMAP_n_doublet_calls.png"),
                      p_n_calls, width = 9, height = 7, dpi = 300)
  }
  
  score_cols <- c("scDblFinder_Score", "scds_Score", "scrublet_Score")
  available_scores <- score_cols[score_cols %in% colnames(obj@meta.data)]
  
  if (length(available_scores) > 0) {
    score_plots <- list()
    for (score_col in available_scores) {
      tool_name <- sub("_Score", "", score_col)
      p_score <- FeaturePlot(obj, features = score_col,
                             reduction = "umap", pt.size = 1, raster = TRUE) +
        scale_color_viridis_c(option = "magma", name = "Score") +
        ggtitle(paste(tool_name, "Score")) +
        theme(legend.position = "right",
              plot.title = element_text(hjust = 0.5, face = "bold"))
      score_plots[[tool_name]] <- p_score
    }
    
    if (length(score_plots) > 0) {
      p_scores_combined <- wrap_plots(score_plots, ncol = 2) +
        plot_annotation(
          title = "Doublet Detection Scores by Tool",
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        )
      save_plot_quietly(file.path(doublet_plots_dir, "UMAP_doublet_scores_comparison.png"),
                        p_scores_combined, width = 14, height = 10, dpi = 300)
    }
  }
  
  if ("doublet_consensus" %in% colnames(obj@meta.data)) {
    doublet_summary <- as.data.frame(obj@meta.data) %>%
      select(contains("Droplet")) %>%
      pivot_longer(cols = everything(), names_to = "Tool", values_to = "Call") %>%
      filter(!is.na(Call)) %>%
      group_by(Tool, Call) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Tool) %>%
      mutate(Percentage = 100 * Count / sum(Count),
             Tool = sub("_Droplet", "", Tool))
    
    p_summary <- ggplot(doublet_summary, aes(x = Tool, y = Percentage, fill = Call)) +
      geom_col(position = "stack") +
      geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
                position = position_stack(vjust = 0.5), size = 3.5, color = "white") +
      scale_fill_manual(values = doublet_colors) +
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Doublet Detection Summary Across Tools",
           x = "Tool", y = "Percentage of Cells", fill = "Classification")
    save_plot_quietly(file.path(doublet_plots_dir, "Barplot_doublet_summary.png"),
                      p_summary, width = 10, height = 6, dpi = 300)
  }
  
  log_status("PLOT", sprintf("Doublet detection plots saved to: %s/", basename(doublet_plots_dir)), "SUCCESS")
  
  # ========== CONTINUE WITH ORIGINAL PLOTS ==========
  log_status("PLOT", "Generating QC plots...", "PROGRESS")
  
  qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  p6 <- FeaturePlot(obj, features = qc_feats, reduction = "umap",
                    order = TRUE, raster = TRUE, ncol = 3) & NoLegend()
  save_plot_quietly(file.path(plots_dir, "UMAP_QC_metrics.png"), p6,
                    width = 15, height = 5, dpi = 300)
  
  markers <- intersect(c("EPCAM", "KRT8", "KRT18", "PTPRC", "CD3D", "CD14"), rownames(obj))
  if (length(markers) > 0) {
    p7 <- FeaturePlot(obj, features = markers, reduction = "umap",
                      order = TRUE, raster = TRUE, ncol = 3,
                      min.cutoff = "q5", max.cutoff = "q95") & NoLegend()
    save_plot_quietly(file.path(plots_dir, "UMAP_key_markers.png"), p7,
                      width = 15, height = 10, dpi = 300)
  }
  
  if ("Monaco_main" %in% colnames(obj@meta.data)) {
    log_status("PLOT", "Generating SingleR annotation plots...", "PROGRESS")
    
    p_monaco <- DimPlot(obj, reduction = "umap", group.by = "Monaco_main",
                        pt.size = 0.2, raster = TRUE, label = TRUE, repel = TRUE) +
      ggtitle("UMAP: Monaco Immune Cell Types")
    save_plot_quietly(file.path(plots_dir, "UMAP_Monaco_main.png"), p_monaco,
                      width = 10, height = 8, dpi = 300)
    
    p_monaco_pruned <- DimPlot(obj, reduction = "umap", group.by = "Monaco_main_pruned",
                               pt.size = 0.2, raster = TRUE, label = TRUE, repel = TRUE) +
      ggtitle("UMAP: Monaco (pruned)")
    save_plot_quietly(file.path(plots_dir, "UMAP_Monaco_main_pruned.png"), p_monaco_pruned,
                      width = 10, height = 8, dpi = 300)
  }
  
  if ("HPCA_main" %in% colnames(obj@meta.data)) {
    p_hpca <- DimPlot(obj, reduction = "umap", group.by = "HPCA_main",
                      pt.size = 0.2, raster = TRUE, label = TRUE, repel = TRUE) +
      ggtitle("UMAP: HPCA Cell Types")
    save_plot_quietly(file.path(plots_dir, "UMAP_HPCA_main.png"), p_hpca,
                      width = 10, height = 8, dpi = 300)
    
    p_hpca_pruned <- DimPlot(obj, reduction = "umap", group.by = "HPCA_main_pruned",
                             pt.size = 0.2, raster = TRUE, label = TRUE, repel = TRUE) +
      ggtitle("UMAP: HPCA (pruned)")
    save_plot_quietly(file.path(plots_dir, "UMAP_HPCA_main_pruned.png"), p_hpca_pruned,
                      width = 10, height = 8, dpi = 300)
  }
  
  v1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "Subject_fixed", pt.size = 0, ncol = 3)
  save_plot_quietly(file.path(plots_dir, "Violin_QC_by_Subject.png"), v1,
                    width = 14, height = 5, dpi = 300)
  
  if ("TimePoint_fixed" %in% colnames(obj@meta.data)) {
    comp_data <- as.data.frame(obj@meta.data) %>%
      mutate(cluster = as.character(cluster)) %>%
      filter(TimePoint_fixed %in% c("Before", "After")) %>%
      group_by(TimePoint_fixed, cluster) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(TimePoint_fixed) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup()
    
    p9 <- ggplot(comp_data, aes(x = cluster, y = frac, fill = TimePoint_fixed)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7) +
      scale_y_continuous(labels = percent_format()) +
      theme_bw(12) +
      labs(title = "Cluster composition: Before vs After",
           x = "Cluster", y = "Fraction", fill = "TimePoint")
    save_plot_quietly(file.path(plots_dir, "Barplot_cluster_composition_BeforeAfter.png"), p9,
                      width = 12, height = 6, dpi = 300)
  }
  
  count_data <- as.data.frame(table(obj$Lane, obj$Subject_fixed))
  colnames(count_data) <- c("Lane", "Subject", "Count")
  
  p10 <- ggplot(count_data, aes(x = Subject, y = Lane, fill = Count)) +
    geom_tile() +
    geom_text(aes(label = scales::comma(Count)), size = 3) +
    scale_fill_gradient(low = "white", high = "#3B82F6") +
    theme_bw(12) +
    labs(title = "Cell counts (Lane × Subject)", x = "Subject", y = "Lane")
  save_plot_quietly(file.path(plots_dir, "Heatmap_counts_Lane_vs_Subject.png"), p10,
                    width = 8.5, height = 5.5, dpi = 300)
  
  log_status("PLOT", "All plots generated successfully", "SUCCESS")
  
  return(obj)
}
################################################################################
# MAIN PIPELINE
################################################################################

main <- function() {
  section_header("SINGLE-CELL RNA-SEQ ANALYSIS PIPELINE")
  
  start_time <- Sys.time()
  log_status("INIT", sprintf("Project: %s", PROJECT), "INFO")
  log_status("INIT", sprintf("Output directory: %s", OUT_DIR), "INFO")
  
  obj_raw <- stage1_load_and_integrate()
  obj_raw <- stage2_add_metadata(obj_raw)
  
  raw_file <- file.path(OUT_DIR, "object_raw_with_metadata.rds")
  log_status("SAVE", "Saving raw object with metadata...", "PROGRESS")
  ##saveRDS(obj_raw, raw_file)
  log_status("SAVE", sprintf("Saved: %s (%.1f MB)", 
                             basename(raw_file),
                             file.size(raw_file) / 1024^2), "SUCCESS")
  
  obj_clean <- stage3_filter_and_harmonize(obj_raw)
  obj_clean <- stage4_dimension_reduction(obj_clean)
  obj_filtered <- stage5_qc_filtering(obj_clean)
  
  log_status("DIMRED", "Re-running dimension reduction on filtered cells...", "PROGRESS")
  obj_filtered <- stage4_dimension_reduction(obj_filtered)
  
  obj_filtered <- stage6_cd45_gating(obj_filtered)
  obj_filtered <- stage7_singler_annotation(obj_filtered)
  obj_filtered <- stage8_cluster_annotation(obj_filtered)
  obj_filtered <- stage9_generate_plots(obj_filtered)
  
  final_file <- file.path(OUT_DIR, "object_final_annotated.rds")
  log_status("SAVE", "Saving final annotated object...", "PROGRESS")
  saveRDS(obj_filtered, final_file)
  log_status("SAVE", sprintf("Saved: %s (%.1f MB)",
                             basename(final_file),
                             file.size(final_file) / 1024^2), "SUCCESS")
  
  section_header("PIPELINE SUMMARY")
  
  
  # Save doublet summary
  if ("doublet_consensus" %in% colnames(obj_filtered@meta.data)) {
    doublet_summary_file <- file.path(OUT_DIR, "doublet_detection_summary.csv")
    
    doublet_summary_df <- obj_filtered@meta.data %>%
      select(contains("Droplet"), contains("_Score"), 
             doublet_consensus, doublet_confidence, 
             n_doublet_calls, n_tools_called) %>%
      mutate(cell_id = rownames(obj_filtered@meta.data))
    
    fwrite(doublet_summary_df, doublet_summary_file)
    log_status("SAVE", sprintf("SAVED: %s", basename(doublet_summary_file)), "SUCCESS")
    
    cat(sprintf("\n"))
    cat(sprintf("[SUMMARY] DOUBLET DETECTION SUMMARY:\n"))
    
    consensus_table <- table(obj_filtered$doublet_consensus)
    for (call_type in names(consensus_table)) {
      pct <- 100 * consensus_table[call_type] / ncol(obj_filtered)
      cat(sprintf("   %s: %s cells (%.1f%%)\n", 
                  call_type, format(consensus_table[call_type], big.mark = ","), pct))
    }
    
    cat(sprintf("\n"))
    cat(sprintf("   Per-tool doublet rates:\n"))
    for (tool_col in c("DoubletDetection_Droplet", "scDblFinder_Droplet", 
                       "scds_Droplet", "scrublet_Droplet")) {
      if (tool_col %in% colnames(obj_filtered@meta.data)) {
        tool_name <- sub("_Droplet", "", tool_col)
        n_dbl <- sum(obj_filtered@meta.data[[tool_col]] == "doublet", na.rm = TRUE)
        n_total <- sum(!is.na(obj_filtered@meta.data[[tool_col]]))
        pct_dbl <- 100 * n_dbl / n_total
        cat(sprintf("     - %s: %.1f%% (%s/%s)\n", 
                    tool_name, pct_dbl, 
                    format(n_dbl, big.mark = ","),
                    format(n_total, big.mark = ",")))
      }
    }
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat(sprintf("\n"))
  cat(sprintf("[COMPLETE] ANALYSIS COMPLETE\n"))
  cat(sprintf("   Total runtime: %.1f minutes\n", as.numeric(elapsed)))
  cat(sprintf("\n"))
  cat(sprintf("[FILES] OUTPUT FILES:\n"))
  cat(sprintf("   1. Raw object: %s\n", basename(raw_file)))
  cat(sprintf("   2. Final object: %s\n", basename(final_file)))
  cat(sprintf("   3. Annotations: cluster_annotations_comprehensive.csv\n"))
  cat(sprintf("   4. Markers: markers/ directory (with caching)\n"))
  cat(sprintf("   5. SingleR: singler_results/\n"))
  cat(sprintf("   6. Plots: plots/ directory\n"))
  cat(sprintf("\n"))
  cat(sprintf("[DATA] DATA SUMMARY:\n"))
  cat(sprintf("   Raw cells: %s\n", format(ncol(obj_raw), big.mark = ",")))
  cat(sprintf("   Filtered cells: %s\n", format(ncol(obj_filtered), big.mark = ",")))
  cat(sprintf("   Clusters identified: %d\n", length(unique(obj_filtered$cluster))))
  cat(sprintf("   Genes analyzed: %s\n", format(nrow(obj_filtered), big.mark = ",")))
  
  if ("CD45pos" %in% colnames(obj_filtered@meta.data)) {
    cd45_pct <- 100 * sum(obj_filtered$CD45pos) / ncol(obj_filtered)
    cat(sprintf("   CD45+ cells: %.1f%%\n", cd45_pct))
  }
  
  if ("TimePoint_fixed" %in% colnames(obj_filtered@meta.data)) {
    tp_table <- table(obj_filtered$TimePoint_fixed)
    cat(sprintf("\n"))
    cat(sprintf("[TIMEPOINT]TIMEPOINT DISTRIBUTION:\n"))
    for (tp in names(tp_table)) {
      pct <- 100 * tp_table[tp] / ncol(obj_filtered)
      cat(sprintf("   %s: %s cells (%.1f%%)\n", 
                  tp, format(tp_table[tp], big.mark = ","), pct))
    }
  }
  
  cat(sprintf("\n"))
  cat(sprintf("[OK] All outputs saved to: %s\n", OUT_DIR))
  cat(sprintf("\n"))
  cat(sprintf("[TIP] TIP: Set ANNOTATION_CONFIG$force_recompute = TRUE to regenerate markers\n"))
  cat(sprintf("\n"))
  
  log_status("COMPLETE", "Pipeline finished successfully!", "SUCCESS")
  
  invisible(obj_filtered)
}

################################################################################
# RUN PIPELINE
################################################################################

if (!interactive()) {
  main()
}

