#!/bin/bash
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --job-name=seg_vis
#SBATCH --output=seg_vis_%A_%a.out
#SBATCH --error=seg_vis_%A_%a.err

echo "Loading conda environment..."
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Define paths
INPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/outputs/segmentation_regions"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/outputs/allele_stacks/all_samples/sample_segmentation_haplotypes"
SCRIPT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/python"

mkdir -p $OUTPUT_DIR


python ${SCRIPT_DIR}/IGV-segmentation.py "$INPUT_DIR" "$OUTPUT_DIR" 
