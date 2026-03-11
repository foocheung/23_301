# Set your project root (you already did this)
PRJ=/data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs

mkdir -p "$PRJ/SNPS8/COMBINED_DEMUX"

for i in 1 2 3 4 5 6; do
  RUNID="23-301-3_GEX${i}"

  DEMUX_FILE="$PRJ/SNPS8/DEMUX/${RUNID}_demux/assignments_refined.tsv.gz"  # or assignments.tsv.gz
##  VIREO_DIR="$PRJ/SNPS8/VIREO_GF/${RUNID}/vireo"
  DD_DIR="$PRJ/SNPS8/DOUBLETD/DoubletDetection/${RUNID}"
  SCDblF_DIR="$PRJ/SNPS8/DOUBLETD/scDblFinder/${RUNID}"
  SCDS_DIR="$PRJ/SNPS8/DOUBLETD/scds/${RUNID}"
  SCRUB_DIR="$PRJ/SNPS8/DOUBLETD/Scrublet/${RUNID}"
  OUT_FILE="$PRJ/SNPS8/COMBINED_DEMUX/${RUNID}_combined.tsv"

  # Sanity checks (skip if any key file/dir is missing)
  if [ ! -f "$DEMUX_FILE" ]; then
    echo "[SKIP] Missing $DEMUX_FILE"
    continue
  fi
##  if [ ! -d "$VIREO_DIR" ]; then
##    echo "[SKIP] Missing $VIREO_DIR"
##    continue
##  fi

  echo "[RUN] $RUNID"
  singularity exec -B "$PRJ":"$PRJ" Demuxafy.sif \
    /opt/Demultiplexing_Doublet_Detecting_Docs/scripts/Combine_Results.R \
    -o "$OUT_FILE" \
    -z "$DEMUX_FILE" \
##    -v "$VIREO_DIR" \
    -t "$DD_DIR" \
    -n "$SCDblF_DIR" \
    -c "$SCDS_DIR" \
    -r "$SCRUB_DIR" \
    -b Demuxalot   \
    -m MajoritySinglet
done

# (Optional) concatenate all per-run combined TSVs into one file
MASTER="$PRJ/SNPS8/COMBINED_DEMUX/23-301-3_ALL_runs_combined.tsv"
# take header from the first file, then append others without header; add a run column
first=1
> "$MASTER"
for i in 1 2 3 4 5 6; do
  RUNID="23-301-3_GEX${i}"
  f="$PRJ/SNPS8/COMBINED_DEMUX/${RUNID}_combined.tsv"
  [ ! -f "$f" ] && continue
  if [ $first -eq 1 ]; then
    awk -v R="$RUNID" 'NR==1{$0=$0"\trun"} NR==1 || NR>1{$0=$0"\t"R}1' OFS='\t' "$f" > "$MASTER"
    first=0
  else
    awk -v R="$RUNID" 'NR>1{$0=$0"\t"R} NR>1' OFS='\t' "$f" >> "$MASTER"
  fi
done

echo "[DONE] Wrote $MASTER"

