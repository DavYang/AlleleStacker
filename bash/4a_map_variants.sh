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
CNV_VCF="${BASE_DIR}/hificnv_postprocessed_vcf/0/SC1.GRCh38.hificnv.ploidy_fixed.vcf.gz"
SV_VCF="${BASE_DIR}/pbsv_jasminesv_postprocessed_vcf/0/SC1.GRCh38.pbsv.jasminesv.ploidy_fixed.vcf.gz"
TR_VCF="${BASE_DIR}/trgt_vcf/0/0/SC1.GRCh38.trgt.adotto_repeats.hg38.vcf.gz"

# Set methylation BED files
H1_BED="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/filtered_consensus_regions/min_sample-1/filtered_consensus_H1.bed"
H2_BED="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/filtered_consensus_regions/min_sample-1/filtered_consensus_H2.bed"

# Function to check and index VCF if needed
check_and_index_vcf() {
    local vcf_file=$1
    local vcf_type=$2
    
    if [[ ! -f "$vcf_file" ]]; then
        echo "Error: VCF file not found: $vcf_file"
        return 1
    fi
    
    # Check if index exists and is newer than VCF
    if [[ ! -f "${vcf_file}.tbi" ]] || [[ "${vcf_file}.tbi" -ot "$vcf_file" ]]; then
        echo "Index not found or outdated for ${vcf_type} VCF. Creating index..."
        
        # Check if file is gzipped
        if [[ "$vcf_file" != *.gz ]]; then
            echo "VCF file is not gzipped. Compressing first..."
            bgzip -c "$vcf_file" > "${vcf_file}.gz"
            vcf_file="${vcf_file}.gz"
        fi
        
        # Create index with progress feedback
        echo "Indexing ${vcf_type} VCF..."
        if tabix -p vcf "$vcf_file" 2>/dev/null; then
            echo "✓ Successfully indexed ${vcf_type} VCF"
        else
            echo "Error: Failed to index ${vcf_type} VCF"
            return 1
        fi
    else
        echo "✓ Valid index found for ${vcf_type} VCF"
    fi
    
    return 0
}

# Create timestamped output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)


# Start logging
exec 1> >(tee "${OUTPUT_DIR}/run.log") 2>&1

echo "Starting variant mapping at $(date)"
echo "Output directory: $OUTPUT_DIR"
echo "Using $(nproc) CPU cores"
echo "Memory limit: 64GB"

# Validate and index VCF files
echo "Checking and indexing VCF files..."
declare -A vcf_files=(
    ["Small variants"]="$SMALL_VCF"
    ["CNV"]="$CNV_VCF"
    ["SV"]="$SV_VCF"
    ["TR"]="$TR_VCF"
)

for vcf_type in "${!vcf_files[@]}"; do
    vcf_path="${vcf_files[$vcf_type]}"
    if ! check_and_index_vcf "$vcf_path" "$vcf_type"; then
        echo "Error processing ${vcf_type} VCF"
        exit 1
    fi
done

# Function to process a single haplotype
process_haplotype() {
    local hap=$1
    local bed_file=$2
    local output_prefix=$3
    
    echo "Processing ${hap} haplotype..."
    
    # Calculate optimal thread count (75% of available CPUs)
    local thread_count=$(( $(nproc) * 3 / 4 ))
    
    $PYTHON/variant_mapper.py \
        --bed "$bed_file" \
        --output "${output_prefix}_${hap}_variants.tsv" \
        --small-vcf "$SMALL_VCF" \
        --cnv-vcf "$CNV_VCF" \
        --sv-vcf "$SV_VCF" \
        --tr-vcf "$TR_VCF" 
        

    if [[ $? -ne 0 ]]; then
        echo "Error: Variant mapping failed for ${hap} haplotype"
        return 1
    fi

    echo "Completed processing ${hap} haplotype"
    return 0
}

# Process both haplotypes
for haplotype in "H1" "H2"; do
    bed_var="${haplotype}_BED"
    process_haplotype "$haplotype" "${!bed_var}" "${OUTPUT_DIR}/variant_mapping"
    if [[ $? -ne 0 ]]; then
        echo "Error processing ${haplotype}"
        exit 1
    fi
done

# Create summary report
echo "Creating summary report..."
{
    echo "Variant Mapping Summary"
    echo "======================"
    echo "Run timestamp: ${TIMESTAMP}"
    echo ""
    echo "System Information:"
    echo "------------------"
    echo "CPU cores: $(nproc)"
    echo "Memory: $(free -h | awk '/Mem:/ {print $2}')"
    echo "Processing time: $SECONDS seconds"
    echo ""
    echo "Input Files:"
    echo "------------"
    for vcf_type in "${!vcf_files[@]}"; do
        echo "${vcf_type} VCF: ${vcf_files[$vcf_type]}"
        echo "  Size: $(du -h "${vcf_files[$vcf_type]}" | cut -f1)"
        echo "  Index: $(du -h "${vcf_files[$vcf_type]}.tbi" 2>/dev/null | cut -f1 || echo "Not found")"
    done
    echo ""
    echo "Results:"
    echo "--------"
    
    # Count variants by type for each haplotype
    for hap in "H1" "H2"; do
        result_file="${OUTPUT_DIR}/variant_mapping_${hap}_variants.tsv"
        if [[ -f "$result_file" ]]; then
            echo "${hap} Results:"
            total_variants=$(($(wc -l < "$result_file") - 1))
            echo "Total variants: ${total_variants}"
            echo "File size: $(du -h "$result_file" | cut -f1)"
            echo "By type:"
            for type in "SMALL" "CNV" "SV" "TR"; do
                count=$(grep -c "^.*\t${type}\t" "$result_file" || true)
                echo "- ${type}: ${count}"
            done
            echo ""
        fi
    done
    
} > "${OUTPUT_DIR}/mapping_summary.txt"

# Compress results
echo "Compressing results..."
tar -czf "${OUTPUT_DIR}.tar.gz" "${OUTPUT_DIR}"

echo "Variant mapping completed at $(date)"
echo "Total runtime: $SECONDS seconds"
echo "Results available in: ${OUTPUT_DIR}"
echo "Compressed archive: ${OUTPUT_DIR}.tar.gz"
echo "Summary report: ${OUTPUT_DIR}/mapping_summary.txt"

# Performance metrics
echo "Performance metrics:"
echo "-------------------"
echo "CPU usage: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}')%"
echo "Memory usage: $(free -h | awk '/Mem:/ {print $3}')"