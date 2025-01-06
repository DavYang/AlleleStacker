#!/bin/bash
# run_methylation_analysis.sh
#SBATCH --job-name=segmentation_regions
#SBATCH --output=logs/segmentation_regions_%A_%a.out
#SBATCH --error=logs/segmentation_regions_%A_%a.err
#SBATCH --array=1-21
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

# Activate conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Get command line arguments
INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Validate input arguments
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: sbatch run_methylation_analysis.sh <input_dir> <output_dir>"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Get the list of sample names excluding non-sample directories
SAMPLES=($(ls -d "$INPUT_DIR"/SPM* | grep -v "log_files\|regions_by_label\|segmentation_scripts" | xargs -n1 basename | sort))

# Debug print
echo "Total samples: ${#SAMPLES[@]}"
echo "Available samples: ${SAMPLES[*]}"

# Validate array task ID
if [ "$SLURM_ARRAY_TASK_ID" -gt "${#SAMPLES[@]}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) is greater than number of samples (${#SAMPLES[@]})"
    exit 1
fi

# Get the current sample based on the array task ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID - 1]}
echo "Processing sample: $SAMPLE"

# Run the Python script for the current sample
python ./python/segmentation_results.py "$INPUT_DIR" "$OUTPUT_DIR" "$SAMPLE"

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully for $SAMPLE. Results are in $OUTPUT_DIR"
else
    echo "An error occurred during the analysis for $SAMPLE."
    exit 1
fi