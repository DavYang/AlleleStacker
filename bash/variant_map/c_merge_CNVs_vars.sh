#!/bin/bash
#SBATCH --job-name=vcf_merge
#SBATCH --time=infinite
#SBATCH --partition=unlimited
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --output=logs/CNV_merge_%j.out
#SBATCH --error=logs/CNV_merge_%j.err

# Set paths
FILEPATHS="$1"
OUT_DIR="$2"

MERGED_VCF="${OUT_DIR}/merged_CNVs_GRCh38.vcf.gz"
INT_VCF="${OUT_DIR}/merged_CNVs_GRCh38.filtered.vcf.gz"
INT_VCF_2="${OUT_DIR}/merged_CNVs_GRCh38.filtered.renamed.vcf.gz"
FINAL_VCF="${OUT_DIR}/merged_CNVs_GRCh38.qc.vcf.gz"

STATS_FILE="${OUT_DIR}/variant_filtering_stats.txt"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"


source /public/apps/conda3/etc/profile.d/conda.sh
conda activate bcftools


# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Start logging
log "Starting VCF processing pipeline"

# Get initial variant count
INITIAL_COUNT=$(bcftools merge \
    --force-samples \
    --merge both \
    --file-list "$FILEPATHS" \
    --output-type z \
    --threads 5 \
    --output "$MERGED_VCF" && bcftools view -H "$MERGED_VCF" | wc -l)

log "Initial variant count: $INITIAL_COUNT"
tabix -p vcf "$MERGED_VCF"

# Filter, rename, and sort variants
bcftools filter -i 'QUAL>=20 && FILTER="PASS"' "$MERGED_VCF" -O z -o "$INT_VCF"
tabix -p vcf "$INT_VCF"

bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT:%SVLEN' "$INT_VCF" -O z -o "$INT_VCF_2"
tabix -p vcf "$INT_VCF_2"

bcftools sort "$INT_VCF_2" -O z -o "$FINAL_VCF"
tabix -p vcf "$FINAL_VCF"

# Get final count
FINAL_COUNT=$(bcftools view -H "$FINAL_VCF" | wc -l)
REMOVED_COUNT=$((INITIAL_COUNT - FINAL_COUNT))

# Write statistics
{
    echo "VCF Processing Summary"
    echo "====================="
    echo "Initial variants: $INITIAL_COUNT"
    echo "Final variants: $FINAL_COUNT"
    echo "Variants removed: $REMOVED_COUNT"
    echo "Percent retained: $(awk "BEGIN {printf \"%.2f%%\", ($FINAL_COUNT/$INITIAL_COUNT)*100}")"
} > "$STATS_FILE"

log "Processing complete. Results in $FINAL_VCF"
cat "$STATS_FILE"