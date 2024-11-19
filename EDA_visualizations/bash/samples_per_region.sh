#!/bin/bash
#SBATCH --job-name=methyl_dist
#SBATCH --output=./log/methyl_dist_%j.out
#SBATCH --error=./log/methyl_dist_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=quick

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Define the input files
H1_INPUT="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/filtered_consensus_regions/min_sample-1/filtered_consensus_H1.bed"
H2_INPUT="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/filtered_consensus_regions/min_sample-1/filtered_consensus_H2.bed"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/EDA_visualizations/outputs/samples_per_region/min_sample-1"

# H1_INPUT="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/filtered_consensus_regions/min_sample-2/filtered_consensus_H1.bed"
# H2_INPUT="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/filtered_consensus_regions/min_sample-2/filtered_consensus_H2.bed"
# OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/EDA_visualizations/outputs/samples_per_region/min_sample-2"

OUTPUT_PREFIX="${OUTPUT_DIR}/samples_per_region"
mkdir -p ${OUTPUT_DIR}


# Run the analysis script
python /gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/EDA_visualizations/python/samples_per_region.py \
    --h1-input "${H1_INPUT}" \
    --h2-input "${H2_INPUT}" \
    --output-prefix "${OUTPUT_PREFIX}"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully"
    echo "Results can be found at:"
    echo "  - Plot: ${OUTPUT_PREFIX}_histograms.png"
    echo "  - Statistics: ${OUTPUT_PREFIX}_statistics.txt"
else
    echo "Error occurred during analysis"
    exit 1
fi