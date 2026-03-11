#!/usr/bin/env Rscript
# Runs Demuxalot on 23-301-3_GEX1..6 using FILTERED barcodes (gz -> tsv)

# -------- USER CONFIG (edit base_dir if different) --------
base_dir <- "/data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs"
sif      <- file.path(base_dir, "SNPS", "Demuxafy.sif")
#vcf      <- file.path(base_dir, "SNPS", "snps.vcf")
#vcf      <- file.path(base_dir, "SNPS", "snps.noMissing.exonic.20x.vcf.gz")

vcf      <- file.path(base_dir, "SNPS", "snps.test.vcf.gz")
donors   <- file.path(base_dir, "SNPS", "ind.txt")     # or "8" if using count
donors_is_file <- TRUE
out_root <- file.path(base_dir, "SNPS9", "DEMUX")
filtered_rel <- "outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
# ---------------------------------------------------------

fail <- function(...) { message(sprintf(...)); quit(status = 1) }
ok   <- function(...) { cat(sprintf(...), "\n") }

# sanity
for (p in list(base_dir, sif, vcf)) if (!file.exists(p) && !dir.exists(p)) fail("Missing: %s", p)
if (donors_is_file && !file.exists(donors)) fail("Donor file missing: %s", donors)

# discover BAMs
bam_paths <- Sys.glob(file.path(base_dir, "23-301-3_GEX?", "outs", "possorted_genome_bam.bam"))
bai_paths <- Sys.glob(file.path(base_dir, "23-301-3_GEX?", "outs", "possorted_genome_bam.bam.bai"))
if (!length(bam_paths)) fail("No BAMs found under %s", file.path(base_dir, "23-301-3_GEX?"))
if (length(bam_paths) != length(bai_paths)) fail("Each BAM must have a .bai index")

samples <- lapply(bam_paths, function(bam) {
  sample_dir  <- dirname(dirname(bam))                 # .../23-301-3_GEXx
  sample_name <- basename(sample_dir)
  bar_gz      <- file.path(sample_dir, filtered_rel)
  if (!file.exists(bar_gz)) fail("Missing filtered barcodes for %s: %s", sample_name, bar_gz)
  list(sample_name = sample_name, bam = bam, bar_gz = bar_gz)
})

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

gz_to_tsv <- function(gz, tsv, sample_name) {
  # Remove any partial file
  if (file.exists(tsv)) file.remove(tsv)
  # Use gzip -dc to stream; requires gzip in PATH on the node
  cmd <- sprintf("gzip -dc '%s' > '%s'", gz, tsv)
  ok("[%s] Decompressing filtered barcodes -> %s", sample_name, tsv)
  status <- system(cmd)
  if (status != 0) fail("[%s] gzip -dc failed (status %d) on %s", sample_name, status, gz)
  if (!file.exists(tsv) || file.info(tsv)$size == 0) fail("[%s] barcodes.tsv empty or missing: %s", sample_name, tsv)
}

for (s in samples) {
  sample_out <- file.path(out_root, paste0(s$sample_name, "_demux"))
  log_dir    <- file.path(sample_out, "logs")
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir,    recursive = TRUE, showWarnings = FALSE)

  bar_tsv   <- file.path(sample_out, "barcodes.tsv")
  log_out   <- file.path(log_dir, "demux_stdout.log")
  log_err   <- file.path(log_dir, "demux_stderr.log")

  # Decompress filtered barcodes in a robust way
  gz_to_tsv(s$bar_gz, bar_tsv, s$sample_name)

  ok("=== Running %s ===", s$sample_name)
  ok("BAM: %s", s$bam); ok("VCF: %s", vcf); ok("Barcodes: %s", bar_tsv); ok("Out: %s", sample_out)

  args <- c("exec", sif, "Demuxalot.py",
            "-a", s$bam,
            "-v", vcf,
            "-b", bar_tsv,
            "-o", sample_out,
            "-r", "True")
  if (donors_is_file) args <- c(args, "-n", donors) else args <- c(args, "-n", as.character(donors))

  ok("singularity %s", paste(args, collapse = " "))

  status <- system2("singularity", args, stdout = log_out, stderr = log_err)
  if (status != 0) {
    message(sprintf("Demuxalot FAILED for %s (status=%d). Logs:\n  %s\n  %s",
                    s$sample_name, status, log_out, log_err))
    next
  }
  n_out <- length(list.files(sample_out, all.files = TRUE))
  ok("Done: %s  (files in out: %d)", s$sample_name, n_out)
}

ok("All samples processed. Root out: %s", out_root)

