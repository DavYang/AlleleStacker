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

# Create directory structure
mkdir -p "${output_dir}/"{regions,summaries}

# Process a single bed file
process_bed() {
    local bed_file="$1"
    local sample_name="$2"
    
    [ ! -f "$bed_file" ] && return
    
    # Process each unique label
    awk -v sample="$sample_name" -v outdir="$regions_dir" -v temp="$temp_dir" '
    NR > 1 && $4 != "summary_label" {
        label = $4
        size = $3 - $2
        
        # Create directory if needed
        if (!dirs[label]++) {
            system("mkdir -p " outdir "/" label)
        }
        
        # Write to label file
        if (!written[label]++) {
            print "chrom\tstart\tend\tsummary_label\tsize" > outdir "/" label "/" sample "_" label ".bed"
        }
        print $1, $2, $3, $4, size >> outdir "/" label "/" sample "_" label ".bed"
        
        # Count regions
        counts[label]++
    }
    END {
        for (label in counts) {
            print label "\t" counts[label] > temp "/" sample "_counts.txt"
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
        process_bed "$bed_file" "${spm_dir}.${hap}"
    done
done

# Generate summary
{
    # Header
    echo -n "Label"
    for spm_dir in SPM*/; do
        [ -d "$spm_dir" ] && echo -n -e "\t${spm_dir%/}.hap1\t${spm_dir%/}.hap2"
    done
    echo
    
    # Get all unique labels and sort them
    find . -name "*.meth_regions.bed" -exec awk 'NR > 1 && $4 != "summary_label" {print $4}' {} \; | 
    sort -u | while read label; do
        echo -n "$label"
        for spm_dir in SPM*/; do
            [ ! -d "$spm_dir" ] && continue
            spm=${spm_dir%/}
            for hap in hap1 hap2; do
                count=$(grep -P "^$label\t" "$temp_dir/${spm}.${hap}_counts.txt" 2>/dev/null | cut -f2)
                echo -n -e "\t${count:-0}"
            done
        done
        echo
    done
} > "${output_dir}/summaries/summary.txt"

rm -r "$temp_dir"
echo "Processing complete. Results in $output_dir"