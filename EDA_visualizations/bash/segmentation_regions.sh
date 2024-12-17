#!/bin/bash
#SBATCH --job-name=segmentation_regions
#SBATCH --output=logs/segmentation_regions_%A_%a.out
#SBATCH --error=logs/segmentation_regions_%A_%a.err
#SBATCH --array=1-21
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

INPUT_DIR="$1"
OUTPUT_DIR="$2"

mkdir -p $OUTPUT_DIR

# Debug exactly where we are and what's here
echo "Current directory:"
pwd
echo -e "\nDirectory contents:"
ls -la
echo -e "\nParent directory contents:"
ls -la ..

# Get the list of sample names excluding non-sample directories
SAMPLES=($(ls -d SPM* | grep -v "log_files\|regions_by_label\|segmentation_scripts" | sort))

# Debug print
echo -e "\nFound samples: ${SAMPLES[@]}"
echo "Total samples: ${#SAMPLES[@]}"
echo "Processing array task: $SLURM_ARRAY_TASK_ID"

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