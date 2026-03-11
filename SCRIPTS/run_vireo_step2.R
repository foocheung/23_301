#!/usr/bin/env Rscript
# Genotype-free Vireo on existing cellSNP outputs (step 3 only).
# Inputs : SNPS3/CELLSNP_ONLY/<sample>/cellsnp/
# Outputs: SNPS3/VIREO_GF/<sample>/vireo/

# -------- USER CONFIG --------
base_dir     <- "/data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs"
sif          <- file.path(base_dir, "SNPS", "Demuxafy.sif")
cellsnp_root <- file.path(base_dir, "SNPS8", "CELLSNP_ONLY")  # matches your ls
out_root     <- file.path(base_dir, "SNPS8", "VIREO_GF")
bind_path    <- base_dir
threads      <- 8
donors_N     <- as.integer(Sys.getenv("VIREO_N", unset = "8"))  # override with env var
rand_seed    <- 12345
# --------------------------------

fail <- function(...) { message(sprintf(...)); quit(status = 1) }
ok   <- function(...) { cat(sprintf(...), "\n") }

# sanity
for (p in list(base_dir, sif, cellsnp_root)) if (!file.exists(p) && !dir.exists(p)) fail("Missing: %s", p)

# discover samples
sample_dirs <- Sys.glob(file.path(cellsnp_root, "23-301-3_GEX?"))
if (!length(sample_dirs)) fail("No cellSNP dirs under %s", cellsnp_root)

# helper to run singularity
singexec <- function(args, stdout = NULL, stderr = NULL) {
  args <- c("exec", "--bind", bind_path, sif, args)
  system2("singularity", args, stdout = stdout, stderr = stderr)
}

# detect support for --randSeed (older vireo may not have it)
has_randSeed <- local({
  out <- tryCatch(system2("singularity",
                          c("exec", "--bind", bind_path, sif, "vireo", "-h"),
                          stdout = TRUE, stderr = TRUE), error = function(e) "")
  any(grepl("--randSeed", out, fixed = TRUE))
})

ok("Vireo K (nDonor) = %d; threads = %d; randSeed supported = %s",
   donors_N, threads, if (has_randSeed) "yes" else "no")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (sd in sample_dirs) {
  sample_name <- basename(sd)
  cs_dir      <- file.path(sd, "cellsnp")
  if (!dir.exists(cs_dir)) { message("Skipping ", sample_name, " (no cellsnp dir)"); next }

  # quick input checks: need AD/DP and base VCF (gz or not)
  req <- c("cellSNP.base.vcf", "cellSNP.base.vcf.gz",
           "cellSNP.tag.AD.mtx", "cellSNP.tag.AD.mtx.gz",
           "cellSNP.tag.DP.mtx", "cellSNP.tag.DP.mtx.gz")
  have <- file.exists(file.path(cs_dir, req))
  if (!any(have[1:2])) fail("[%s] Missing cellSNP.base.vcf(.gz)", sample_name)
  if (!any(have[3:4])) fail("[%s] Missing cellSNP.tag.AD.mtx(.gz)", sample_name)
  if (!any(have[5:6])) fail("[%s] Missing cellSNP.tag.DP.mtx(.gz)", sample_name)

  sample_out <- file.path(out_root, sample_name, "vireo")
  log_dir    <- file.path(out_root, sample_name, "logs")
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir,    recursive = TRUE, showWarnings = FALSE)

  vi_stdout <- file.path(log_dir, "vireo_stdout.log")
  vi_stderr <- file.path(log_dir, "vireo_stderr.log")

  ok("=== Vireo (genotype-free) :: %s ===", sample_name)
  ok("cellSNP dir: %s", cs_dir); ok("Out: %s", sample_out)

  v_args <- c("vireo",
              "-c", cs_dir,
              "-o", sample_out,
              "-N", as.character(donors_N),
              "-p", as.character(threads))
  if (has_randSeed) v_args <- c(v_args, "--randSeed", as.character(rand_seed))

  status <- singexec(v_args, stdout = vi_stdout, stderr = vi_stderr)
  if (status != 0) {
    message(sprintf("Vireo FAILED for %s (exit %d). Logs:\n  %s\n  %s",
                    sample_name, status, vi_stdout, vi_stderr))
    next
  }

  assign_tsv <- if (file.exists(file.path(sample_out, "donor_ids.tsv.gz")))
                  file.path(sample_out, "donor_ids.tsv.gz") else
                  file.path(sample_out, "donor_ids.tsv")

  ok("Done: %s  (assignments: %s%s)",
     sample_name, basename(assign_tsv),
     if (file.exists(assign_tsv)) "" else " [not found]")
}

ok("All samples processed with genotype-free Vireo. Root out: %s", out_root)

