#!/usr/bin/env Rscript
# Dropulation (per-sample): annotate BAM -> assign -> doublets -> final calls
# Inputs expected (per sample): outs/possorted_genome_bam.bam, outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# Also needs: donor VCF (+ tabix index), donor list (INDS), GTF
# Outputs: SNPS3/DROPULATION/<sample>/{assignments.tsv.gz, likelihoods.tsv.gz, updated_assignments.tsv.gz, ...}

## ---------------- USER CONFIG ----------------
base_dir  <- "/data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs"
sif       <- file.path(base_dir, "SNPS", "Demuxafy.sif")
vcf       <- file.path(base_dir, "SNPS", "test.vcf.gz")
inds      <- file.path(base_dir, "SNPS", "ind.txt")      # donor IDs present in VCF (one per line)
gtf       <- file.path(base_dir, "SNPS", "genes.gtf")                        # <-- set this (ideally the same GTF used by Cell Ranger)
cell_tag  <- "CB"                                        # SAM tag for cell barcode
umi_tag   <- "UB"                                        # SAM tag for UMI
threads   <- 12
bind_path <- base_dir
out_root  <- file.path(base_dir, "SNPS9", "DROPULATION")
filtered_rel <- "outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
## ---------------------------------------------

fail <- function(...) { message(sprintf(...)); quit(status = 1) }
ok   <- function(...) { cat(sprintf(...), "\n") }

for (p in list(base_dir, sif, vcf, inds, gtf)) if (!file.exists(p) && !dir.exists(p)) fail("Missing: %s", p)

# find samples
bam_paths <- Sys.glob(file.path(base_dir, "23-301-3_GEX?", "outs", "possorted_genome_bam.bam"))
if (!length(bam_paths)) fail("No BAMs found under %s", file.path(base_dir, "23-301-3_GEX?"))
if (length(bam_paths) != length(Sys.glob(paste0(bam_paths, ".bai")))) fail("Each BAM must have a .bai index")

# helpers
singexec <- function(args, log_out, log_err) {
  args <- c("exec", "--bind", bind_path, sif, args)
  system2("singularity", args, stdout = log_out, stderr = log_err)
}
ensure_tabix <- function(vcfgz, log_dir) {
  if (!file.exists(paste0(vcfgz, ".tbi")) && !file.exists(paste0(vcfgz, ".csi"))) {
    ok("[*] Indexing VCF: %s", vcfgz)
    st <- singexec(c("tabix", "-p", "vcf", vcfgz),
                   file.path(log_dir, "tabix_stdout.log"), file.path(log_dir, "tabix_stderr.log"))
    if (st != 0) fail("tabix failed on %s", vcfgz)
  }
}
gz_to_tsv <- function(gz, tsv, sample_name) {
  if (file.exists(tsv)) file.remove(tsv)
  cmd <- sprintf("gzip -dc '%s' > '%s'", gz, tsv)
  ok("[%s] Writing barcodes.tsv", sample_name); st <- system(cmd)
  if (st != 0 || !file.exists(tsv) || file.info(tsv)$size == 0) fail("[%s] barcodes.tsv failed", sample_name)
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
ensure_tabix(vcf, out_root)

for (bam in bam_paths) {
  sample_dir  <- dirname(dirname(bam))
  sample_name <- basename(sample_dir)
  bar_gz      <- file.path(sample_dir, filtered_rel)
  if (!file.exists(bar_gz)) fail("[%s] Missing %s", sample_name, bar_gz)

  sample_out <- file.path(out_root, sample_name)
  log_dir    <- file.path(sample_out, "logs")
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir,    recursive = TRUE, showWarnings = FALSE)

  bar_tsv <- file.path(sample_out, "barcodes.tsv")
  gz_to_tsv(bar_gz, bar_tsv, sample_name)

  tagged_bam <- file.path(sample_out, "possorted_genome_bam_dropulation_tag.bam")

  # 1) Annotate BAM with gene function (if not already present)
  if (!file.exists(tagged_bam)) {
    ok("[%s] TagReadWithGeneFunction -> %s", sample_name, basename(tagged_bam))
    st <- singexec(c("TagReadWithGeneFunction",
                     paste0("ANNOTATIONS_FILE=", gtf),
                     paste0("INPUT=", bam),
                     paste0("OUTPUT=", tagged_bam)),
                   file.path(log_dir, "tag_stdout.log"), file.path(log_dir, "tag_stderr.log"))
    if (st != 0) { message(sprintf("[%s] Tagging failed (see logs).", sample_name)); next }
  } else ok("[%s] Reusing existing %s", sample_name, basename(tagged_bam))

  # 2) Assign (per-cell singlet donor likelihood)
  assign_tsv <- file.path(sample_out, "assignments.tsv.gz")
  ok("[%s] Dropulation_AssignCellsToSamples.py", sample_name)
  st <- singexec(c("Dropulation_AssignCellsToSamples.py",
                   "--CELL_BC_FILE", bar_tsv,
                   "--INPUT_BAM",   tagged_bam,
                   "--OUTPUT",      assign_tsv,
                   "--VCF",         vcf,
                   "--SAMPLE_FILE", inds,
                   "--CELL_BARCODE_TAG", cell_tag,
                   "--MOLECULAR_BARCODE_TAG", umi_tag,
                   "--VCF_OUTPUT",  file.path(sample_out, "assignment.vcf"),
                   "--MAX_ERROR_RATE", "0.05"),
                 file.path(log_dir, "assign_stdout.log"),
                 file.path(log_dir, "assign_stderr.log"))
  if (st != 0) { message(sprintf("[%s] Assignment failed.", sample_name)); next }

  # 3) Doublet likelihoods
  like_tsv <- file.path(sample_out, "likelihoods.tsv.gz")
  ok("[%s] DetectDoublets", sample_name)
  st <- singexec(c("DetectDoublets",
                   "--CELL_BC_FILE", bar_tsv,
                   "--INPUT_BAM",    tagged_bam,
                   "--OUTPUT",       like_tsv,
                   "--VCF",          vcf,
                   "--SINGLE_DONOR_LIKELIHOOD_FILE", assign_tsv,
                   "--SAMPLE_FILE",  inds,
                   "--CELL_BARCODE_TAG", cell_tag,
                   "--MOLECULAR_BARCODE_TAG", umi_tag,
                   "--MAX_ERROR_RATE", "0.05"),
                 file.path(log_dir, "doublet_stdout.log"),
                 file.path(log_dir, "doublet_stderr.log"))
  if (st != 0) { message(sprintf("[%s] DetectDoublets failed.", sample_name)); next }

  # 4) Final calls
  final_tsv <- file.path(sample_out, "updated_assignments.tsv.gz")
  ok("[%s] dropulation_call.R", sample_name)
  st <- singexec(c("dropulation_call.R",
                   "--assign", assign_tsv,
                   "--doublet", like_tsv,
                   "--out", final_tsv),
                 file.path(log_dir, "call_stdout.log"),
                 file.path(log_dir, "call_stderr.log"))
  if (st != 0) { message(sprintf("[%s] dropulation_call failed.", sample_name)); next }

  ok("Done: %s  -> %s", sample_name, final_tsv)
}
ok("All samples processed. Root out: %s", out_root)

