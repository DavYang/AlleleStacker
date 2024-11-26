#!/bin/bash
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --job-name=meth_vis
#SBATCH --output=meth_vis_%A_%a.out
#SBATCH --error=meth_vis_%A_%a.err
#SBATCH --array=0-1
#SBATCH --mail-type=END
#SBATCH --mail-user=david.yang@einsteinmed.edu

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Function for error handling
error_handler() {
    echo "Error occurred in script at line: ${1}"
    echo "Exit code: ${2}"
}
trap 'error_handler ${LINENO} $?' ERR

# Load conda environment
echo "Loading conda environment..."
source /public/apps/conda3/etc/profile.d/conda.sh
conda activate anc_vig

# Define paths
INPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/segmentation_regions"
OUTPUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/allele_stacks/all_samples"
SCRIPT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/python"
SAMPLE_LIST="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/sample_lists/all_samples.txt"

mkdir -p $OUTPUT_DIR

# Array of haplotypes
HAPLOTYPES=("H1" "H2")
HAPLOTYPE=${HAPLOTYPES[$SLURM_ARRAY_TASK_ID]}

echo "Starting processing at $(date)"
echo "Processing haplotype: $HAPLOTYPE"

# Run the processing script for this haplotype
python ${SCRIPT_DIR}/IGV-all-samples.py \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --haplotype "$HAPLOTYPE" \
    --sample_list "$SAMPLE_LIST" \
    --threads 20

echo "Python processing completed at $(date)"

# Set up output file names
BED_FILE="${OUTPUT_DIR}/igv_cohort_${HAPLOTYPE,,}.bed"

if [[ -f "$BED_FILE" ]]; then
    echo "Processing output files..."
    
    # First sort the BED file (excluding track lines)
    echo "Sorting BED data..."
    # Save track lines
    grep "^track" "$BED_FILE" > "${BED_FILE}.tracks"
    # Sort data lines (excluding track lines)
    grep -v "^track" "$BED_FILE" | \
        sort -k1,1 -k2,2n > "${BED_FILE}.sorted"
    
    # Combine track lines and sorted data
    echo "Combining track definitions and sorted data..."
    cat "${BED_FILE}.tracks" "${BED_FILE}.sorted" > "${BED_FILE}.temp"
    mv "${BED_FILE}.temp" "$BED_FILE"
    
    # Clean up temporary files
    rm "${BED_FILE}.tracks" "${BED_FILE}.sorted"
    
    # Compress the file
    echo "Compressing BED file..."
    bgzip -f "$BED_FILE"
    
    # Index the compressed file
    echo "Creating tabix index..."
    tabix -p bed -f "${BED_FILE}.gz"
    
    echo "Output files created:"
    echo "  ${BED_FILE}.gz (compressed BED)"
    echo "  ${BED_FILE}.gz.tbi (index)"
    
    # Verify files exist
    if [[ -f "${BED_FILE}.gz" && -f "${BED_FILE}.gz.tbi" ]]; then
        echo "Successfully created compressed and indexed files"
    else
        echo "Error: Failed to create output files"
        exit 1
    fi
else
    echo "Error: BED file not found: $BED_FILE"
    exit 1
fi

echo "Job completed successfully at $(date)"
