#!/bin/bash
set -euo pipefail

# -------------------------------
# CONFIGURATION
# -------------------------------
# Base paths
BASE_DIR="/data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs"
SNPS_DIR="${BASE_DIR}/SNPS9/DOUBLETD"

# Software container
SIF="./SNPS/Demuxafy.sif"

# Create log directory
LOGDIR="${SNPS_DIR}/logs"
mkdir -p "$LOGDIR"

# Number of expected doublets (adjust per sample if needed)
N_DOUB=3200

# Array of all samples
SAMPLES=("23-301-3_GEX1" "23-301-3_GEX2" "23-301-3_GEX3" "23-301-3_GEX4" "23-301-3_GEX5" "23-301-3_GEX6")

# -------------------------------
# RUNNING scDblFinder FOR ALL SAMPLES (PARALLEL)
# -------------------------------
echo "========================================"
echo "BATCH: Running scDblFinder for all samples (parallel)"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    (
        echo "Processing ${SAMPLE} with scDblFinder..."
        COUNTS="${BASE_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix"
        SCDBLFINDER_OUTDIR="${SNPS_DIR}/scDblFinder/${SAMPLE}"
        mkdir -p "$SCDBLFINDER_OUTDIR"
        
        if [ ! -d "$COUNTS" ]; then
            echo "✗ ERROR: Matrix not found for ${SAMPLE}: $COUNTS"
            exit 1
        fi
        
        singularity exec "$SIF" scDblFinder.R \
          -o "$SCDBLFINDER_OUTDIR" \
          -t "$COUNTS" \
          > "${LOGDIR}/${SAMPLE}_scDblFinder.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ scDblFinder completed for ${SAMPLE}"
        else
            echo "✗ scDblFinder failed for ${SAMPLE} - check ${LOGDIR}/${SAMPLE}_scDblFinder.log"
        fi
    ) &
done

wait
echo "scDblFinder batch completed!"
echo ""

# -------------------------------
# RUNNING DoubletDetection FOR ALL SAMPLES (PARALLEL)
# -------------------------------
echo "========================================"
echo "BATCH: Running DoubletDetection for all samples (parallel)"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    (
        echo "Processing ${SAMPLE} with DoubletDetection..."
        COUNTS="${BASE_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix"
        DOUBLETDETECTION_OUTDIR="${SNPS_DIR}/DoubletDetection/${SAMPLE}"
        mkdir -p "$DOUBLETDETECTION_OUTDIR"
        
        if [ ! -d "$COUNTS" ]; then
            echo "✗ ERROR: Matrix not found for ${SAMPLE}"
            exit 1
        fi
        
        singularity exec "$SIF" DoubletDetection.py \
          -m "$COUNTS" \
          -o "$DOUBLETDETECTION_OUTDIR" \
          > "${LOGDIR}/${SAMPLE}_DoubletDetection.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ DoubletDetection completed for ${SAMPLE}"
        else
            echo "✗ DoubletDetection failed for ${SAMPLE} - check ${LOGDIR}/${SAMPLE}_DoubletDetection.log"
        fi
    ) &
done

wait
echo "DoubletDetection batch completed!"
echo ""

# -------------------------------
# RUNNING scds FOR ALL SAMPLES (PARALLEL)
# -------------------------------
echo "========================================"
echo "BATCH: Running scds for all samples (parallel)"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    (
        echo "Processing ${SAMPLE} with scds..."
        COUNTS="${BASE_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix"
        SCDS_OUTDIR="${SNPS_DIR}/scds/${SAMPLE}"
        mkdir -p "$SCDS_OUTDIR"
        
        if [ ! -d "$COUNTS" ]; then
            echo "✗ ERROR: Matrix not found for ${SAMPLE}"
            exit 1
        fi
        
        singularity exec "$SIF" scds.R \
          -o "$SCDS_OUTDIR" \
          -t "$COUNTS" \
          > "${LOGDIR}/${SAMPLE}_scds.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ scds completed for ${SAMPLE}"
        else
            echo "✗ scds failed for ${SAMPLE} - check ${LOGDIR}/${SAMPLE}_scds.log"
        fi
    ) &
done

wait
echo "scds batch completed!"
echo ""

# -------------------------------
# RUNNING Scrublet FOR ALL SAMPLES (PARALLEL)
# -------------------------------
echo "========================================"
echo "BATCH: Running Scrublet for all samples (parallel)"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    (
        echo "Processing ${SAMPLE} with Scrublet..."
        COUNTS="${BASE_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix"
        SCRUBLET_OUTDIR="${SNPS_DIR}/Scrublet/${SAMPLE}"
        mkdir -p "$SCRUBLET_OUTDIR"
        
        if [ ! -d "$COUNTS" ]; then
            echo "✗ ERROR: Matrix not found for ${SAMPLE}"
            exit 1
        fi
        
        singularity exec "$SIF" Scrublet.py \
          -m "$COUNTS" \
          -o "$SCRUBLET_OUTDIR" \
          > "${LOGDIR}/${SAMPLE}_Scrublet.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ Scrublet completed for ${SAMPLE}"
        else
            echo "✗ Scrublet failed for ${SAMPLE} - check ${LOGDIR}/${SAMPLE}_Scrublet.log"
        fi
    ) &
done

wait
echo "Scrublet batch completed!"
echo ""

# -------------------------------
# RUNNING DoubletDecon FOR ALL SAMPLES (PARALLEL)
# -------------------------------
echo "========================================"
echo "BATCH: Running DoubletDecon for all samples (parallel)"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    (
        echo "Processing ${SAMPLE} with DoubletDecon..."
        
        # DoubletDecon requires a Seurat object
        SEURAT_RDS="${BASE_DIR}/${SAMPLE}/seurat_object.rds"
        DOUBLETDECON_OUTDIR="${SNPS_DIR}/DoubletDecon/${SAMPLE}"
        mkdir -p "$DOUBLETDECON_OUTDIR"
        
        if [ ! -f "$SEURAT_RDS" ]; then
            echo "⚠ WARNING: Seurat object not found for ${SAMPLE}, skipping DoubletDecon"
            exit 0
        fi
        
        singularity exec "$SIF" DoubletDecon.R \
          -o "$DOUBLETDECON_OUTDIR" \
          -s "$SEURAT_RDS" \
          > "${LOGDIR}/${SAMPLE}_DoubletDecon.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ DoubletDecon completed for ${SAMPLE}"
        else
            echo "✗ DoubletDecon failed for ${SAMPLE} - check ${LOGDIR}/${SAMPLE}_DoubletDecon.log"
        fi
    ) &
done

wait
echo "DoubletDecon batch completed!"
echo ""

# -------------------------------
# RUNNING DoubletFinder FOR ALL SAMPLES (PARALLEL)
# -------------------------------
echo "========================================"
echo "BATCH: Running DoubletFinder for all samples (parallel)"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    (
        echo "Processing ${SAMPLE} with DoubletFinder..."
        
        # DoubletFinder requires a Seurat object
        SEURAT_RDS="${BASE_DIR}/${SAMPLE}/seurat_object.rds"
        DOUBLETFINDER_OUTDIR="${SNPS_DIR}/DoubletFinder/${SAMPLE}"
        mkdir -p "$DOUBLETFINDER_OUTDIR"
        
        if [ ! -f "$SEURAT_RDS" ]; then
            echo "⚠ WARNING: Seurat object not found for ${SAMPLE}, skipping DoubletFinder"
            exit 0
        fi
        
        singularity exec "$SIF" DoubletFinder.R \
          -o "$DOUBLETFINDER_OUTDIR" \
          -s "$SEURAT_RDS" \
          -c TRUE \
          -d "$N_DOUB" \
          > "${LOGDIR}/${SAMPLE}_DoubletFinder.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ DoubletFinder completed for ${SAMPLE}"
        else
            echo "✗ DoubletFinder failed for ${SAMPLE} - check ${LOGDIR}/${SAMPLE}_DoubletFinder.log"
        fi
    ) &
done

wait
echo "DoubletFinder batch completed!"
echo ""

# -------------------------------
# RUNNING Solo FOR ALL SAMPLES (PARALLEL)
# -------------------------------
echo "========================================"
echo "BATCH: Running Solo for all samples (parallel)"
echo "========================================"

for SAMPLE in "${SAMPLES[@]}"; do
    (
        echo "Processing ${SAMPLE} with Solo..."
        COUNTS="${BASE_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix"
        SOLO_OUTDIR="${SNPS_DIR}/Solo/${SAMPLE}"
        mkdir -p "$SOLO_OUTDIR"
        
        # Check if JSON exists for this sample
        JSON="${BASE_DIR}/${SAMPLE}/samplesheet.json"
        if [ ! -f "$JSON" ]; then
            echo "⚠ WARNING: JSON not found for ${SAMPLE}, skipping Solo"
            exit 0
        fi
        
        if [ ! -d "$COUNTS" ]; then
            echo "✗ ERROR: Matrix not found for ${SAMPLE}"
            exit 1
        fi
        
        singularity exec "$SIF" solo \
          -o "$SOLO_OUTDIR" \
          -e "$N_DOUB" \
          -j "$JSON" \
          -d "$COUNTS" \
          > "${LOGDIR}/${SAMPLE}_Solo.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "✓ Solo completed for ${SAMPLE}"
        else
            echo "✗ Solo failed for ${SAMPLE} - check ${LOGDIR}/${SAMPLE}_Solo.log"
        fi
    ) &
done

wait
echo "Solo batch completed!"
echo ""

# -------------------------------
# SUMMARY
# -------------------------------
echo "========================================"
echo "ALL BATCHES COMPLETED!"
echo "========================================"
echo "Logs available in: ${LOGDIR}"
echo ""
echo "Summary by tool:"
echo "- scDblFinder: ${SNPS_DIR}/scDblFinder/"
echo "- DoubletDetection: ${SNPS_DIR}/DoubletDetection/"
echo "- scds: ${SNPS_DIR}/scds/"
echo "- Scrublet: ${SNPS_DIR}/Scrublet/"
echo "- DoubletDecon: ${SNPS_DIR}/DoubletDecon/"
echo "- DoubletFinder: ${SNPS_DIR}/DoubletFinder/"
echo "- Solo: ${SNPS_DIR}/Solo/"
echo ""
echo "Check individual logs for any failures:"
ls -lh "${LOGDIR}"
echo "========================================"
