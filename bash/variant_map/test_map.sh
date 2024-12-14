#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --array=1,2

# Usage function
usage() {
    echo "Usage: $0 -b BASE_DIR -o OUTPUT_DIR [-h]"
    echo "  -b: Base directory containing AlleleStacker outputs"
    echo "  -o: Output directory for results"
    echo "  -h: Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "b:o:h" opt; do
    case $opt in
        b) BASE_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [[ -z "$BASE_DIR" || -z "$OUTPUT_DIR" ]]; then
    usage
fi

# Exit on error
set -e

# Map array task to haplotype
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]]; then
    HAP="H1"
else
    HAP="H2"
fi

# Define paths based on base directory
VARIANT_DIR="${BASE_DIR}/outputs/variant_mapping/merged_variants"
PYTHON_DIR="${BASE_DIR}/python/variant_map"
BED_FILE="${BASE_DIR}/outputs/candidate_regions/filtered/filtered_candidate_${HAP}.bed"


# Create output directory
mkdir -p "$OUTPUT_DIR"

# Load conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Start logging
exec 1> >(tee "${OUTPUT_DIR}/${HAP}.log") 2>&1

echo "Starting variant mapping for ${HAP} at $(date)"
echo "Using $(nproc) CPU cores"

# Check input files
if [[ ! -f "$BED_FILE" ]]; then
    echo "Error: BED file not found: $BED_FILE"
    exit 1
fi

# Run variant mapping
python3 "${PYTHON_DIR}/test_map2.py" \
    --bed "$BED_FILE" \
    --haplotype "$HAP" \
    --output-prefix "${OUTPUT_DIR}/${HAP}" \
    --small-vcf "${VARIANT_DIR}/small_variants/merged_phased_small_GRCh38.deepvariant.qc.vcf.gz" \
    --cnv-vcf "${VARIANT_DIR}/copy_number_variants/merged_CNVs_GRCh38.qc.vcf.gz" \
    --sv-vcf "${VARIANT_DIR}/structural_variants/merged_phased_SVs_GRCh38.qc.vcf.gz" \
    --threads ${SLURM_CPUS_PER_TASK}

echo "Completed variant mapping for ${HAP} at $(date)"

# Create simple summary
echo -e "\nResults Summary:" > "${OUTPUT_DIR}/${HAP}_summary.txt"
echo "Total variants: $(( $(wc -l < "${OUTPUT_DIR}/${HAP}_variants.tsv") - 1 ))" >> "${OUTPUT_DIR}/${HAP}_summary.txt"
echo "Runtime: $SECONDS seconds" >> "${OUTPUT_DIR}/${HAP}_summary.txt"

echo "All done! Results in ${OUTPUT_DIR}"