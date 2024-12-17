#!/bin/bash
#SBATCH --job-name=regions-per-chr
#SBATCH --output=logs/regions_%A_%a.out
#SBATCH --error=logs/regions_%A_%a.err
#SBATCH --partition=quick
#SBATCH --mem=16G
#SBATCH --time=02:00:00

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_prefix>"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR=$1
OUTPUT_DIR=$2


mkdir -p "$OUTPUT_DIR"

# Run the Python script with the provided arguments
python ./python/regions_by_chrom.py "$INPUT_DIR" "$OUTPUT_DIR/"