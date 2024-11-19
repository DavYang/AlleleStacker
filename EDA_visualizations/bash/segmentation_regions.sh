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

INPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/segmentation_regions"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/EDA_visualizations/outputs/segmentation_EDA/region_counts"

# Get the list of sample names
SAMPLES=($(ls $INPUT_DIR/*.meth_regions.bed | xargs -n 1 basename | sed 's/\.meth_regions\.bed//'))

# Get the current sample based on the array task ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID - 1]}

# Run the Python script for the current sample
python ./python/segmentation_results.py "$INPUT_DIR" "$OUTPUT_DIR" "$SAMPLE"

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully for $SAMPLE. Results are in $OUTPUT_DIR"
else
    echo "An error occurred during the analysis for $SAMPLE."
    exit 1
fi
