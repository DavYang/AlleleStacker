#!/bin/bash
#SBATCH --job-name=region_plots
#SBATCH --output=logs/region_plots_%A.out
#SBATCH --error=logs/region_plots_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_dir output_dir"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2

mkdir -p "$OUTPUT_DIR"
mkdir -p logs


python python/segmentation/plot_region_distribution_counts.py --base_dir "$INPUT_DIR" --output_dir "$OUTPUT_DIR"
