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
# This assumes all samples have files in H1_M
SAMPLES=($(ls "$1"/regions/H1_M/SPM*.bed | xargs -n1 basename | sed 's/_H1_M.bed//g' | sort -u))
CURRENT_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$CURRENT_INDEX]}

# Run the Python script with the correct paths
python ./python/segmentation/plot_unmethylated_distribution.py \
    --input-dir "$1" \
    --output-dir "$2" \
    --sample-name "$SAMPLE"
