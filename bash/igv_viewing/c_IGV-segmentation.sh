#!/bin/bash
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=4:00:00
#SBATCH --job-name=seg_vis
#SBATCH --output=logs/seg_vis_%A_%a.out
#SBATCH --error=logs/seg_vis_%A_%a.err

echo "Loading conda environment..."
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Define paths
INPUT_DIR="$1"
OUTPUT_DIR="$2"
SCRIPT_DIR="./python/igv_viewing"

mkdir -p $OUTPUT_DIR


python ${SCRIPT_DIR}/IGV-segmentation.py "$INPUT_DIR" "$OUTPUT_DIR" 
