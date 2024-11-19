#!/bin/bash
#SBATCH --job-name=check
#SBATCH --output=logs/check_%A_%a.out
#SBATCH --error=logs/check_%A_%a.err
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick


python ./python/check_consensus.py \
    --consensus-file /gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/filtered_consensus_regions/filtered_consensus_H1_filtered.bed \
    --regions-dir /gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/segmentation_regions/regions_by_label \
    --haplotype H1 \
    --output-file test_region_analysis.txt