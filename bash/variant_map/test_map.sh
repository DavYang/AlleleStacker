#!/bin/bash
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --job-name=test
#SBATCH --output=var_test_%j.out
#SBATCH --error=var_test_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --array=1,2

# Exit on error
set -e

# Load conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Map array task to haplotype
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]]; then
    HAP="H1"
    BED_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-19-24/outputs/filtered_consensus_regions/min_sample-1/filtered_consensus_H1.bed"
else
    HAP="H2"
    BED_FILE="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-19-24/outputs/filtered_consensus_regions/min_sample-1/filtered_consensus_H2.bed"
fi

# Set paths
VARIANT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/merged_variants"
PYTHON_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/python"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/variant_mapping"

# Create timestamped output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="${OUTPUT_DIR}/test_run_${TIMESTAMP}"
mkdir -p "$RUN_DIR"

# Start logging
exec 1> >(tee "${RUN_DIR}/test_${HAP}.log") 2>&1

echo "Starting test variant mapping for ${HAP} at $(date)"
echo "Using ${BED_FILE}"

# Run quick test
python3 "${PYTHON_DIR}/test_map2.py" \
    --bed "$BED_FILE" \
    --haplotype "$HAP" \
    --output "${RUN_DIR}/test_${HAP}_variants.tsv" \
    --small-vcf "${VARIANT_DIR}/small_variants/merged_phased_small_GRCh38.deepvariant.qc.vcf.gz" \
    --cnv-vcf "${VARIANT_DIR}/copy_number_variants/merged_CNVs_GRCh38.qc.vcf.gz" \
    --sv-vcf "${VARIANT_DIR}/strutural_variants/merged_phased_SVs_GRCh38.qc.vcf.gz" \
    --sample-size 50
    # --tr-vcf "${VARIANT_DIR}/trgt_repeat_vcf/merged_phased_repeats_GRCh38.vcf.gz"

echo "Test complete for ${HAP} at $(date)"
echo "Results in ${RUN_DIR}/test_${HAP}_variants.tsv"

# Generate quick stats
echo -e "\nQuick statistics:"
wc -l "${RUN_DIR}/test_${HAP}_variants.tsv"
echo "Variant counts by type:"
tail -n +2 "${RUN_DIR}/test_${HAP}_variants.tsv" | cut -f5 | sort | uniq -c