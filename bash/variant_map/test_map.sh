#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=2-00:00:00
#SBATCH --job-name=var_map
#SBATCH --output=logs/var_map_%A_%a.out
#SBATCH --error=logs/var_map_%A_%a.err
#SBATCH --cpus-per-task=20
#SBATCH --array=1,2

# Exit on error
set -e

# Load conda environment
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Get input parameters
if [[ ${SLURM_ARRAY_TASK_ID} == 1 ]]; then
    HAP="H1"
else
    HAP="H2"
fi

# Get input directories from command line arguments
if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <bed_dir> <variant_dir> <python_dir> <output_dir>"
    exit 1
fi

BED_DIR=$1
VARIANT_DIR=$2
PYTHON_DIR=$3
OUTPUT_DIR=$4
BED_FILE="${BED_DIR}/filtered_candidate_${HAP}.bed"

# Create timestamped output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="${OUTPUT_DIR}/run_${TIMESTAMP}"
mkdir -p "$RUN_DIR"
mkdir -p "logs"

# Start logging
exec 1> >(tee "${RUN_DIR}/${HAP}.log") 2>&1

echo "Starting variant mapping for ${HAP} at $(date)"
echo "Output directory: $RUN_DIR"
echo "Using $(nproc) CPU cores"
echo "Memory limit: 64GB"

# Check input files
if [[ ! -f "$BED_FILE" ]]; then
    echo "Error: BED file not found: $BED_FILE"
    exit 1
fi

# Run variant mapping
echo "Running variant mapping for ${HAP}..."
python3 "${PYTHON_DIR}/variant_mapper-3.py" \
    --bed "$BED_FILE" \
    --haplotype "$HAP" \
    --output-prefix "${RUN_DIR}/${HAP}" \
    --small-vcf "${VARIANT_DIR}/small_variants/merged_phased_small_GRCh38.deepvariant.qc.vcf.gz" \
    --cnv-vcf "${VARIANT_DIR}/copy_number_variants/merged_CNVs_GRCh38.qc.vcf.gz" \
    --sv-vcf "${VARIANT_DIR}/structural_variants/merged_phased_SVs_GRCh38.qc.vcf.gz" \
    --tr-vcf "${VARIANT_DIR}/tandem_repeats/merged_TRs_GRCh38.qc.vcf.gz" \
    --threads ${SLURM_CPUS_PER_TASK}

# Check exit status
if [[ $? -ne 0 ]]; then
    echo "Error: Variant mapping failed"
    exit 1
fi

echo "Completed variant mapping for ${HAP} at $(date)"
