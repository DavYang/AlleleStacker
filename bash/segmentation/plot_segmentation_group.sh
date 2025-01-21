#!/bin/bash
#SBATCH --job-name=group_analysis
#SBATCH --output=logs/group_analysis_%j.out
#SBATCH --error=logs/group_analysis_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

mkdir -p "$2"/segmentation_plots_group
PYTHON="./python/segmentation"

# Generate summary table
python $PYTHON/segmentation_group_summary.py "$1" "$2"/segmentation_plots_group/all_samples_summary.csv

# Then create plots
python $PYTHON/segmentation_plots_group.py "$1" "$2"/segmentation_plots_group