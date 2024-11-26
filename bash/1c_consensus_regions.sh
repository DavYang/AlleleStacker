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
INPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/segmentation_regions/regions_by_label"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/consensus_regions/min_sample-1"
OUTPUT_PREFIX="consensus_regions"

MIN_SAMPLES=1
MAX_GAP=500
PYTHON="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/python"
SAMPLE_LIST_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/all_samples.txt"



### Test

# # 3 sample
# OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/consensus_regions/min_sample-1/test_3"
# SAMPLE_LIST_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/3_samples.txt"

# # 6 sample
# OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/consensus_regions/min_sample-1/test_6"
# SAMPLE_LIST_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/6_samples.txt"

# # 9 sample
# OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/consensus_regions/min_sample-1/test_9"
# SAMPLE_LIST_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/9_samples.txt"

# # 12 sample
# OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/consensus_regions/min_sample-1/test_12"
# SAMPLE_LIST_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/12_samples.txt"

# # 15 sample
# OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/consensus_regions/min_sample-1/test_15"
# SAMPLE_LIST_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/15_samples.txt"

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

