#!/bin/bash
#SBATCH --job-name=meth_dist
#SBATCH --output=log/meth_dist_%A_%a.out
#SBATCH --error=log/meth_dist_%A_%a.err
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.yang2@einsteinmed.edu

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig


filtered_CONSENSUS_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/filtered_consensus_regions/min_sample-1"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/EDA_visualizations/outputs/region_seg_calls/min_sample-1"

REGIONS_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/segmentation_regions/regions_by_label"


# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

python python/interactive_region_meth_dist.py \
    --h1 $filtered_CONSENSUS_DIR/filtered_consensus_H1.bed \
    --h2 $filtered_CONSENSUS_DIR/filtered_consensus_H2.bed \
    --output $OUTPUT_DIR \
    --n-regions 50
