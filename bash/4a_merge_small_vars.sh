#!/bin/bash
#SBATCH --job-name=vcf_merge
#SBATCH --time=infinite
#SBATCH --partition=unlimited
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=david.yang2@einsteinmed.edu
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --output=vcf_merge_%j.out
#SBATCH --error=vcf_merge_%j.err

# Set paths
OUT_DIR="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/variant_prep/small_variants"
MERGED_VCF="${OUT_DIR}/merged_phased_small_GRCh38.deepvariant.vcf.gz"
FILEPATHS="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker_11-20-24/outputs/variant_prep/file_paths/filepaths_small.txt"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"


source /public/apps/conda3/etc/profile.d/conda.sh
conda activate bcftools

# Merge VCFs using bcftools
bcftools merge \
    --force-samples \
    --merge both \
    --file-list "$FILEPATHS" \
    --output-type z \
    --threads 20 \
    --output "$MERGED_VCF"

# Index the merged VCF
tabix -p vcf "$MERGED_VCF"

# Print completion message
echo "VCF merge complete. Output: $MERGED_VCF"