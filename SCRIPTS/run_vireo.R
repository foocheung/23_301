#!/usr/bin/env Rscript
# Runs ONLY the first step of the pipeline for 23-301-3_GEX1..6:
# barcodes.gz -> barcodes.tsv -> cellSNP-lite pileup
# (No donor-VCF subsetting and NO Vireo.)
#
# Output: /.../SNPS3/CELLSNP_ONLY/<sample>/cellsnp/

# -------- USER CONFIG (edit base_dir if different) --------
base_dir <- "/data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs"
sif      <- file.path(base_dir, "SNPS", "Demuxafy.sif")
vcf      <- file.path(base_dir, "SNPS", "snps.test.vcf.gz")

out_root <- file.path(base_dir, "SNPS9", "CELLSNP_ONLY")
filtered_rel <- "outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

threads <- 8                 # threads for cellSNP-lite
bind_path <- base_dir        # path to bind into the container (parent of all inputs/outputs)
# ---------------------------------------------------------

fail <- function(...) { message(sprintf(...)); quit(status = 1) }
ok   <- function(...) { cat(sprintf(...), "\n") }

# sanity
for (p in list(base_dir, sif, vcf)) if (!file.exists(p) && !dir.exists(p)) fail("Missing: %s", p)

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
  if (file.exists(tsv)) file.remove(tsv)
  cmd <- sprintf("gzip -dc '%s' > '%s'", gz, tsv)
  ok("[%s] Decompressing filtered barcodes -> %s", sample_name, tsv)
  status <- system(cmd)
  if (status != 0) fail("[%s] gzip -dc failed (status %d) on %s", sample_name, status, gz)
  if (!file.exists(tsv) || file.info(tsv)$size == 0) fail("[%s] barcodes.tsv empty or missing: %s", sample_name, tsv)
}

singexec <- function(args, log_out, log_err) {
  args <- c("exec", "--bind", bind_path, sif, args)
  system2("singularity", args, stdout = log_out, stderr = log_err)
}

ensure_vcf_index <- function(vcfgz, log_dir) {
  tbi <- paste0(vcfgz, ".tbi"); csi <- paste0(vcfgz, ".csi")
  if (!file.exists(tbi) && !file.exists(csi)) {
    ok("[*] Indexing VCF (for cellSNP -R): %s", vcfgz)
    status <- singexec(c("tabix", "-p", "vcf", vcfgz),
                       file.path(log_dir, "tabix_stdout.log"),
                       file.path(log_dir, "tabix_stderr.log"))
    if (status != 0) fail("tabix failed on %s (exit %d)", vcfgz, status)
  }
}

for (s in samples) {
  sample_out <- file.path(out_root, s$sample_name)
  log_dir    <- file.path(sample_out, "logs")
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir,    recursive = TRUE, showWarnings = FALSE)

  bar_tsv   <- file.path(sample_out, "barcodes.tsv")
  log_out   <- file.path(log_dir, "stdout.log")
  log_err   <- file.path(log_dir, "stderr.log")

  # Decompress filtered barcodes
  gz_to_tsv(s$bar_gz, bar_tsv, s$sample_name)

  # Ensure donor VCF is indexed (cellSNP uses -R <vcf>)
  ensure_vcf_index(vcf, log_dir)

  ok("=== Running cellSNP-lite for %s ===", s$sample_name)
  ok("BAM: %s", s$bam); ok("VCF: %s", vcf); ok("Barcodes: %s", bar_tsv); ok("Out: %s", sample_out)

  # 1) cellSNP-lite pileup (mode 1a: given SNP list)
  cellsnp_out <- file.path(sample_out, "cellsnp")
  dir.create(cellsnp_out, showWarnings = FALSE)
  cs_stdout <- file.path(log_dir, "cellsnp_stdout.log")
  cs_stderr <- file.path(log_dir, "cellsnp_stderr.log")

  cs_args <- c("cellsnp-lite",
               "-s", s$bam,
               "-b", bar_tsv,
               "-O", cellsnp_out,
               "-R", vcf,
               "-p", as.character(threads),
               "--minMAF", "0.1",
               "--minCOUNT", "20",
               "--gzip")
  ok("[%s] singularity exec ... cellsnp-lite (out=%s)", s$sample_name, cellsnp_out)
  status <- singexec(cs_args, cs_stdout, cs_stderr)
  if (status != 0) {
    message(sprintf("cellSNP-lite FAILED for %s (status=%d). Logs:\n  %s\n  %s",
                    s$sample_name, status, cs_stdout, cs_stderr))
    next
  }

  # Stop here—no donor VCF subsetting and no Vireo
  ok("Done (first step only): %s", s$sample_name)
}

ok("All samples processed (first step only). Root out: %s", out_root)

