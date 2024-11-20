#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=2-00:00:00
#SBATCH --job-name=var_map
#SBATCH --output=var_map_%A_%a.out
#SBATCH --error=var_map_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=4

# Exit on error
set -e

# Load conda environment
echo "Loading conda environment..."
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Set base paths
BASE_DIR="/gs/gsfs0/shared-lab/greally-lab/David/6_base-seq_SC1/CoLoRSdb/outputs/20240424_102735_colors_main/out"
PYTHON="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/python"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/variant_mapping"

mkdir -p "$OUTPUT_DIR"
# Set input paths
SMALL_VCF="${BASE_DIR}/deepvariant_glnexus_postprocessed_vcf/0/SC1.GRCh38.deepvariant.glnexus.norm.ploidy_fixed.vcf.gz"

# Set methylation BED files
H1_BED="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/filtered_consensus_regions/min_sample-1/filtered_consensus_H1.bed"
H2_BED="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/filtered_consensus_regions/min_sample-1/filtered_consensus_H2.bed"

python $PYTHON/small_variant_mapper.py \
  --bed $H1_BED \
  --vcf $SMALL_VCF \
  --output $OUTPUT_DIR/results.tsv