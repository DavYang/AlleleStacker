#!/bin/bash
#SBATCH --job-name=parse
#SBATCH --time=2-00:00:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --output=log-parse.%A_%a.out

input_dir="$1"  
output_dir="$2"
regions_dir="${output_dir}/regions"
temp_dir=$(mktemp -d)

# Create directory structure for combined haplotype-methylation states
mkdir -p "${output_dir}/"{regions/{H1_M,H1_U,H2_M,H2_U},summaries}

# Process a single bed file
process_bed() {
    local bed_file="$1"
    local sample_name="$2"
    local hap="$3"  # hap1 or hap2
    
    [ ! -f "$bed_file" ] && return
    
    # Convert hap1/hap2 to H1/H2
    local hap_prefix="${hap/hap/H}"
    
    # Process and relabel regions
    awk -v sample="$sample_name" -v hap="$hap_prefix" -v outdir="$regions_dir" -v temp="$temp_dir" '
    NR > 1 && $4 != "summary_label" {
        # Create combined label (H1_M or H1_U or H2_M or H2_U)
        combined_label = hap "_" $4
        size = $3 - $2
        
        # Write to appropriate file with modified summary_label
        if (!written[combined_label]++) {
            print "chrom\tstart\tend\tsummary_label\tsize" > outdir "/" combined_label "/" sample "_" combined_label ".bed"
        }
        
        # Write data with modified summary_label
        print $1, $2, $3, combined_label, size >> outdir "/" combined_label "/" sample "_" combined_label ".bed"
        
        # Count regions
        counts[combined_label]++
    }
    END {
        # Write counts for all possible states, even if zero
        states[1] = hap "_M"
        states[2] = hap "_U"
        for (i in states) {
            print sample "\t" hap "\t" states[i] "\t" (counts[states[i]] ? counts[states[i]] : 0) >> temp "/" sample "_counts.txt"
        }
    }' OFS='\t' "$bed_file"
}

# Process each sample directory
cd "$input_dir"
for spm_dir in SPM*; do
    [ ! -d "$spm_dir" ] && continue
    echo "Processing $spm_dir..."
    
    # Process both haplotypes
    for hap in hap1 hap2; do
        bed_file="$spm_dir/${spm_dir}.${hap}.meth_regions.bed"
        process_bed "$bed_file" "${spm_dir}" "$hap"
    done
done

# Generate summary with combined haplotype-methylation states
{
    # Header
    echo -e "Sample\tH1_M\tH1_U\tH2_M\tH2_U"
    
    # Process each sample
    for spm_dir in SPM*/; do
        [ ! -d "$spm_dir" ] && continue
        spm=${spm_dir%/}
        echo -n "$spm"
        
        # Create associative array for counts
        declare -A state_counts=([H1_M]=0 [H1_U]=0 [H2_M]=0 [H2_U]=0)
        
        # Read counts file if it exists
        if [ -f "$temp_dir/${spm}_counts.txt" ]; then
            while IFS=$'\t' read -r sample hap state count; do
                state_counts[$state]=$count
            done < "$temp_dir/${spm}_counts.txt"
        fi
        
        # Output counts in fixed order
        for state in H1_M H1_U H2_M H2_U; do
            echo -n -e "\t${state_counts[$state]}"
        done
        echo
    done
} > "${output_dir}/summaries/summary.txt"

rm -r "$temp_dir"
echo "Processing complete. Results in $output_dir"