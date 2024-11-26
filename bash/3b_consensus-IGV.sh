#!/bin/bash
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --job-name=meth_vis
#SBATCH --output=meth_vis_%A_%a.out
#SBATCH --error=meth_vis_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=david.yang@einsteinmed.edu


# Load conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Define directories
INPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/filtered_consensus_regions/min_sample-1"
OUTPUT_DIR=gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs//allele_stacks/filtered_consensus"
SCRIPT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/python"

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
