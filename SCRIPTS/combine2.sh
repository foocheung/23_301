#!/bin/bash
#SBATCH -J combine_all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o combine_all_%j.out
#SBATCH -e combine_all_%j.err

PRJ=/data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs
OUTDIR="$PRJ/SNPS8/COMBINED_DEMUX"
mkdir -p "$OUTDIR"

echo "[START] Combining results for 23-301-3_GEX1-6"

########## Scenario A: Demuxalot + Doublet tools (no Vireo) ##########
for i in 1 2 3 4 5 6; do
  RUN=23-301-3_GEX${i}
  echo "[RUN] $RUN  (Demux-only + doublet tools)"
  singularity exec -B "$PRJ":"$PRJ" Demuxafy.sif \
    /opt/Demultiplexing_Doublet_Detecting_Docs/scripts/Combine_Results.R \
    -o "$OUTDIR/${RUN}_combined_DemuxOnly.tsv" \
    -z "$PRJ/SNPS8/DEMUX/${RUN}_demux/assignments_refined.tsv.gz" \
    -t "$PRJ/SNPS8/DOUBLETD/DoubletDetection/${RUN}" \
    -n "$PRJ/SNPS8/DOUBLETD/scDblFinder/${RUN}" \
    -c "$PRJ/SNPS8/DOUBLETD/scds/${RUN}" \
    -r "$PRJ/SNPS8/DOUBLETD/Scrublet/${RUN}" \
    -b Demuxalot \
    -m MajoritySinglet
done

########## Scenario B: Demuxalot + Vireo only ##########
for i in 1 2 3 4 5 6; do
  RUN=23-301-3_GEX${i}
  echo "[RUN] $RUN  (Demux + Vireo only)"
  singularity exec -B "$PRJ":"$PRJ" Demuxafy.sif \
    /opt/Demultiplexing_Doublet_Detecting_Docs/scripts/Combine_Results.R \
    -o "$OUTDIR/${RUN}_combined_Demux_Vireo.tsv" \
    -z "$PRJ/SNPS8/DEMUX/${RUN}_demux/assignments_refined.tsv.gz" \
    -v "$PRJ/SNPS8/VIREO_GF/${RUN}/vireo" \
    -b Demuxalot \
    -m MajoritySinglet
done

########## Merge helper function ##########
merge_with_run () {
  pattern="$1"; out="$2"
  first=1; : > "$out"
  for i in 1 2 3 4 5 6; do
    RUN=23-301-3_GEX${i}
    f=$(printf "$pattern" "$RUN")
    [ -f "$f" ] || continue
    if [ $first -eq 1 ]; then
      awk -v R="$RUN" 'NR==1{$0=$0"\trun"} {print $0"\t"R}' OFS="\t" "$f" > "$out"
      first=0
    else
      awk -v R="$RUN" 'NR>1{print $0"\t"R}' OFS="\t" "$f" >> "$out"
    fi
  done
  echo "[DONE] $out"
}

########## Merge per-scenario outputs ##########
merge_with_run "$OUTDIR/%s_combined_DemuxOnly.tsv"   "$OUTDIR/23-301-3_ALL_combined_DemuxOnly.tsv"
merge_with_run "$OUTDIR/%s_combined_Demux_Vireo.tsv" "$OUTDIR/23-301-3_ALL_combined_Demux_Vireo.tsv"

########## Quick singlet count summary ##########
echo -e "\n[SUMMARY] Singlet counts per file:"
for f in $OUTDIR/*_combined_*.tsv; do
  c=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){h[$i]=i}} NR>1 && ( ($h["Final_droplet_type"]=="singlet") || ($h["droplet_type"]=="singlet") ){n++} END{print n+0}' "$f")
  printf "%s\t%s\n" "$(basename "$f")" "$c"
done

echo "[ALL DONE]"

