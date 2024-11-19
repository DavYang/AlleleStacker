#!/bin/bash
#SBATCH --job-name=gaps
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=quick
#SBATCH --output=gap.out

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig


HAP1_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/segmentation_regions/regions_by_label/H1_U"
HAP2_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/segmentation_regions/regions_by_label/H2_U"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/EDA_visualizations/segmentation_EDA/region_gaps"
PYTHON="/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/EDA_visualizations/python"

# Run the Python script
python $PYTHON/region_gap.py "$HAP1_DIR" "$HAP2_DIR" "$OUTPUT_DIR"

# Check if the Python script executed successfully
if [ $? -eq 0 ]; then
    echo "Gap analysis completed successfully. Results are in $OUTPUT_DIR"
else
    echo "An error occurred during gap analysis."
    exit 1
fi