#!/bin/bash
#SBATCH --job-name=vcf_merge
#SBATCH --time=infinite
#SBATCH --partition=unlimited
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --output=log.vcf_merge_%j.out
#SBATCH --error=log.vcf_merge_%j.err

# Set paths
FILEPATHS="$1"
OUT_DIR="$2"

MERGED_VCF="${OUT_DIR}/merged_repeats_GRCh38.vcf.gz"
FINAL_VCF="${OUT_DIR}/merged_repeats_GRCh38.renamed.vcf.gz"
STATS_FILE="${OUT_DIR}/variant_merge_stats.txt"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

source /public/apps/conda3/etc/profile.d/conda.sh
conda activate bcftools

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Start logging
log "Starting VCF merging pipeline"

# Merge VCFs and get initial count
INITIAL_COUNT=$(bcftools merge \
    --force-samples \
    --merge both \
    --file-list "$FILEPATHS" \
    --output-type z \
    --threads 5 \
    --output "$MERGED_VCF" && bcftools view -H "$MERGED_VCF" | wc -l)

log "Total variant count: $INITIAL_COUNT"
tabix -p vcf "$MERGED_VCF"

# Rename variants with chromosome position
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' "$MERGED_VCF" -O z -o "$FINAL_VCF"
tabix -p vcf "$FINAL_VCF"

# Write statistics
{
    echo "VCF Merging Summary"
    echo "==================="
    echo "Total variants: $INITIAL_COUNT"
} > "$STATS_FILE"

log "Merging complete. Results in $FINAL_VCF"
cat "$STATS_FILE"
