#!/bin/bash
#SBATCH --job-name=sample_analysis
#SBATCH --output=logs/sample_analysis_%A_%a.out
#SBATCH --error=logs/sample_analysis_%A_%a.err
#SBATCH --array=1-21
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

mkdir -p "$2"/segmentation_plots_sample
mkdir -p logs

SAMPLES=($(ls -d "$1"/SPM* | grep -v "log_files\|regions_by_label\|segmentation_scripts" | xargs -n1 basename | sort))
CURRENT_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$CURRENT_INDEX]}

python python/segmentation_plots_sample.py "$1" "$2"/segmentation_plots_sample "$SAMPLE"