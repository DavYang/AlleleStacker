#!/bin/bash
#SBATCH --job-name=meth_dist
#SBATCH --output=logs/meth_dist_%A_%a.out
#SBATCH --error=logs/meth_dist_%A_%a.err
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick


source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig


INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

python python/interactive_region_meth_dist.py \
    --h1 $INPUT_DIR/filtered_candidate_H1.bed \
    --h2 $INPUT_DIR/filtered_candidate_H2.bed \
    --output $OUTPUT_DIR \
    --n-regions 50
