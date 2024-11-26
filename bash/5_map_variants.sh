#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=2-00:00:00
#SBATCH --job-name=var_map
#SBATCH --output=var_map_%A_%a.out
#SBATCH --error=var_map_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=david.yang2@einsteinmed.edu
#SBATCH --cpus-per-task=20
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
RUN_DIR="${OUTPUT_DIR}/run_${TIMESTAMP}"
mkdir -p "$RUN_DIR"

# Start logging
exec 1> >(tee "${RUN_DIR}/${HAP}.log") 2>&1

echo "Starting variant mapping for ${HAP} at $(date)"
echo "Output directory: $RUN_DIR"
echo "Using $(nproc) CPU cores"
echo "Memory limit: 64GB"


if [[ ! -f "$BED_FILE" ]]; then
    echo "Error: BED file not found: $BED_FILE"
    exit 1
fi

# Run variant mapping
echo "Running variant mapping for ${HAP}..."
python3 "${PYTHON_DIR}/test_map2.py" \
    --bed "$BED_FILE" \
    --haplotype "$HAP" \
    --output-prefix "${RUN_DIR}/${HAP}" \
    --small-vcf "${VARIANT_DIR}/small_variants/merged_phased_small_GRCh38.deepvariant.qc.vcf.gz" \
    --cnv-vcf "${VARIANT_DIR}/copy_number_variants/merged_CNVs_GRCh38.qc.vcf.gz" \
    --sv-vcf "${VARIANT_DIR}/structural_variants/merged_phased_SVs_GRCh38.qc.vcf.gz" \
     --threads ${SLURM_CPUS_PER_TASK}

# Check exit status
if [[ $? -ne 0 ]]; then
    echo "Error: Variant mapping failed"
    exit 1
fi

echo "Completed variant mapping for ${HAP} at $(date)"

# Create summary report
echo "Generating summary report..."
{
    echo "Variant Mapping Run Summary"
    echo "=========================="
    echo "Run timestamp: ${TIMESTAMP}"
    echo "Haplotype: ${HAP}"
    echo ""
    echo "Results:"
    echo "  Output file: ${RUN_DIR}/${HAP}_variants.tsv"
    echo "  Total variants: $(( $(wc -l < "${RUN_DIR}/${HAP}_variants.tsv") - 1 ))"
    echo ""
    echo "Variant counts by type:"
    echo "----------------------"
    for type in small cnv sv tr; do
        count=$(grep -c "^.*\t${type}\t" "${RUN_DIR}/${HAP}_variants.tsv" || true)
        printf "%-6s: %8d\n" "$type" "$count"
    done
    echo ""
    echo "Resource Usage:"
    echo "  Runtime: $SECONDS seconds"
    echo "  Memory peak: $(free -h | awk '/Mem:/ {print $3}')"
    echo "  CPU usage: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}')%"
} > "${RUN_DIR}/${HAP}_summary.txt"

echo "Summary report written to: ${RUN_DIR}/${HAP}_summary.txt"

# Compress results
echo "Compressing results..."
tar -czf "${RUN_DIR}/${HAP}_results.tar.gz" \
    "${RUN_DIR}/${HAP}_variants.tsv" \
    "${RUN_DIR}/${HAP}_summary.txt" \
    "${RUN_DIR}/${HAP}.log"

echo "Results compressed to: ${RUN_DIR}/${HAP}_results.tar.gz"
echo "All done!"