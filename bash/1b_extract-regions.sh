#!/bin/bash
#SBATCH --job-name=parse
#SBATCH --time=00:10:00
#SBATCH --partition=quick
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --output=log-parse.%A_%a.out

# Directory containing the BED files
bed_dir="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/segmentation_regions"

# Base output directory for extracted regions
output_base_dir="/gs/gsfs0/shared-lab/greally-lab/David/AlleleStacker_tests/AlleleStacker/outputs/segmentation_regions/regions_by_label"

# Create the base output directory if it doesn't exist
mkdir -p "$output_base_dir"

# Temporary directory for storing intermediate label counts
temp_dir=$(mktemp -d)

# Function to extract regions based on the label
extract_regions_by_label() {
    local input_bed="$1"
    local output_dir="$2"
    # Get the sample name from the input file
    local sample_name=$(basename "$input_bed" .meth_regions.bed)
    # Get unique labels, excluding "summary_label"
    local labels=$(awk 'NR > 1 && $4 != "summary_label" {print $4}' "$input_bed" | sort | uniq)

    # Extract regions for each label and count the regions
    for label in $labels; do
        # Create the label directory if it doesn't exist
        local label_dir="$output_dir/$label"
        mkdir -p "$label_dir"
        # Output BED file for the current label
        local output_bed="$label_dir/${sample_name}_${label}.bed"
        # Add header and extract regions with the current label, calculating size
        awk -v label="$label" 'BEGIN { print "chrom\tstart\tend\tsummary_label\tsize" } NR > 1 && $4 == label { print $1, $2, $3, $4, $3 - $2 }' OFS='\t' "$input_bed" > "$output_bed"
        # Count the number of regions for the current label
        local count=$(awk -v label="$label" 'NR > 1 && $4 == label { count++ } END { print count }' "$input_bed")
        echo -e "$label\t$count" >> "$temp_dir/${sample_name}_counts.txt"
        echo "Extracting $label regions from $sample_name to $output_bed"
    done
}

# Cycle through each BED file in the directory
for input_bed in "$bed_dir"/*.meth_regions.bed; do
    if [ -f "$input_bed" ]; then
        # Extract regions for each unique label and write to the respective directory
        extract_regions_by_label "$input_bed" "$output_base_dir"
    else
        echo "Input file $input_bed not found."
    fi
done

# Generate the summary file with samples as columns and labels as rows
summary_file="$output_base_dir/summary.txt"
echo -n -e "Label" > "$summary_file"
for input_bed in "$bed_dir"/*.meth_regions.bed; do
    if [ -f "$input_bed" ]; then
        sample_name=$(basename "$input_bed" .meth_regions.bed)
        echo -n -e "\t$sample_name" >> "$summary_file"
    fi
done
echo "" >> "$summary_file"

labels=$(awk 'NR > 1 && $4 != "summary_label" {print $4}' "$bed_dir"/*.meth_regions.bed | sort | uniq)
for label in $labels; do
    echo -n -e "$label" >> "$summary_file"
    for input_bed in "$bed_dir"/*.meth_regions.bed; do
        if [ -f "$input_bed" ]; then
            sample_name=$(basename "$input_bed" .meth_regions.bed)
            count=$(grep -P "^$label\t" "$temp_dir/${sample_name}_counts.txt" | cut -f2)
            echo -n -e "\t$count" >> "$summary_file"
        fi
    done
    echo "" >> "$summary_file"
done

# Clean up temporary directory
rm -r "$temp_dir"

echo "Summary file generated at $summary_file"
