#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --output=consensus_regions.%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=david.yang2@einsteinmed.edu

# Load conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Set your parameters
INPUT_DIR="$1"
OUTPUT_DIR="$2"
OUTPUT_PREFIX="consensus_regions"

MIN_SAMPLES=1
MAX_GAP=500
PYTHON="./python"
SAMPLE_LIST_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/all_samples.txt"



# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Call the Python script
python $PYTHON/consensus_regions.py \
  --input_dir $INPUT_DIR \
  --output_prefix $OUTPUT_PREFIX \
  --output_dir $OUTPUT_DIR \
  --sample-list-file $SAMPLE_LIST_FILE \
  --max-gap $MAX_GAP \
  --min-samples $MIN_SAMPLES

