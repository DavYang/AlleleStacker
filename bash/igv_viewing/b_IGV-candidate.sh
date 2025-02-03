#!/bin/bash
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --job-name=meth_vis
#SBATCH --output=logs/meth_vis_%A_%a.out
#SBATCH --error=logs/meth_vis_%A_%a.err



# Load conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Define directories
INPUT_DIR="$1"
OUTPUT_DIR="$2"
SCRIPT_DIR="python/igv_viewing"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run the Python script
python3 $SCRIPT_DIR/IGV-consensus.py "$INPUT_DIR" "$OUTPUT_DIR"

# Check if the processing was successful
if [ $? -eq 0 ]; then
    echo "Successfully processed BED files"
    echo "Output files are in: $OUTPUT_DIR"
else
    echo "Error processing BED files" >&2
    exit 1
fi
