#!/usr/bin/env python3

import sys
import os
import argparse
from pathlib import Path

def process_bed_file(input_file, output_file):
    """
    Process a BED file to extract chrom, start, end columns and add GREEN color
    """
    print(f"Processing file: {input_file}")
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # Skip header
        header = fin.readline()
        
        # Process each line
        line_count = 0
        for line in fin:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = fields[1]
            end = fields[2]
            
            # Write the modified line with GREEN color
            output_line = f"{chrom}\t{start}\t{end}\tconsensus_region\t0\t.\t{start}\t{end}\t0,255,0\n"
            fout.write(output_line)
            line_count += 1
        
        print(f"Processed {line_count} lines")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process BED files for IGV viewing with GREEN coloring')
    parser.add_argument('input_dir', help='Directory containing input BED files')
    parser.add_argument('output_dir', help='Directory for output BED files')
    
    args = parser.parse_args()
    
    # Convert to Path objects for easier handling
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Specifically look for H1 and H2 files only
    input_files = ['filtered_consensus_H1.bed', 'filtered_consensus_H2.bed']
    
    for filename in input_files:
        input_file = input_dir / filename
        if not input_file.exists():
            print(f"Warning: {input_file} not found")
            continue
            
        # Create output filename
        output_file = output_dir / f"{input_file.stem}_igv{input_file.suffix}"
        
        try:
            process_bed_file(input_file, output_file)
            print(f"Successfully processed {input_file} -> {output_file}")
        except Exception as e:
            print(f"Error processing {input_file}: {str(e)}", file=sys.stderr)

if __name__ == "__main__":
    main()