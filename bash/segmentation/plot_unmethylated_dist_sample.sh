#!/bin/bash
#SBATCH --job-name=unmeth_dist
#SBATCH --output=logs/unmeth_dist_%A_%a.out
#SBATCH --error=logs/unmeth_dist_%A_%a.err
#SBATCH --array=1-21
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Create output and log directories
mkdir -p "$2"
mkdir -p logs

# Get list of unique sample names from the bed files in H1_M directory
# Using absolute paths and proper directory structure
cd "$1/regions/H1_M" || exit 1
SAMPLES=($(ls SPM*.bed | sed 's/_H1_M.bed//g' | sort -u))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "No sample files found in $1/regions/H1_M"
    exit 1
fi

CURRENT_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$CURRENT_INDEX]}

if [ -z "$SAMPLE" ]; then
    echo "No sample found for index $CURRENT_INDEX"
    exit 1
fi

echo "Processing sample: $SAMPLE"

# Run the Python script with the correct paths
# Assuming the Python script is in the project root's python/segmentation directory
PROJECT_ROOT=$(dirname $(dirname $(dirname "$1")))
python "$PROJECT_ROOT/python/segmentation/plot_unmethylated_distribution.py" \
    --input-dir "$1/regions" \
    --output-dir "$2" \
    --sample-name "$SAMPLE"
