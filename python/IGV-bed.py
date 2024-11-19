import os
import sys
import glob
import logging
import pandas as pd
import pybedtools
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def methylation_status_to_rgb(name):
    """Convert methylation status to RGB color"""
    if '_M' in name:
        return "255,0,0"  # Red for methylated
    elif '_U' in name:
        return "0,0,255"  # Blue for unmethylated
    return "0,0,0"  # Black for unknown

def stagger_overlaps(df):
    """Stagger overlapping regions within a sample"""
    df = df.sort_values(by=['chrom', 'start', 'end'])
    current_end = -1
    offset = 0
    for index, row in df.iterrows():
        if row['start'] < current_end:
            offset += 1
        else:
            offset = 0
        df.at[index, 'start'] += offset
        df.at[index, 'end'] += offset
        current_end = row['end']
    return df

def read_sample_list(sample_file):
    """Read sample IDs from a text file"""
    try:
        with open(sample_file, 'r') as f:
            samples = [line.strip() for line in f if line.strip()]
        logging.info(f"Loaded {len(samples)} samples from {sample_file}")
        return samples  # Return as list to maintain order
    except Exception as e:
        logging.error(f"Error reading sample file {sample_file}: {e}")
        sys.exit(1)

def get_sample_id(filename):
    """Extract sample ID from the filename"""
    basename = os.path.basename(filename)
    return basename.replace('.meth_regions.bed', '')

def create_track_lines(sample_id, haplotype):
    """Create IGV track line for a sample"""
    return f'track name="{sample_id}_{haplotype}" description="{sample_id} {haplotype} Methylation" visibility=full useScore=0 itemRgb="On"'

def aggregate_haplotype_regions(input_dir: str, output_dir: str, haplotype: str, sample_list: list = None):
    """
    Aggregate regions and create IGV-compatible BED file with methylation colors and ordered samples.
    """
    logging.info(f"Processing {haplotype} regions...")
    
    # Output file path
    output_file = os.path.join(output_dir, f"igv_methylation_{haplotype.lower()}.bed")
    
    # Get all meth_regions.bed files
    bed_files = glob.glob(os.path.join(input_dir, "*meth_regions.bed"))
    if not bed_files:
        logging.error(f"No *meth_regions.bed files found in {input_dir}")
        return

    # If no sample list provided, create one from available files
    if sample_list is None:
        sample_list = sorted([get_sample_id(f) for f in bed_files])
    
    # Process each sample in order and write to file
    with open(output_file, 'w') as outfile:
        for sample_id in sample_list:
            bed_file = os.path.join(input_dir, f"{sample_id}.meth_regions.bed")
            if not os.path.exists(bed_file):
                logging.warning(f"File not found for sample {sample_id}, skipping")
                continue
                
            try:
                # Write track line for this sample
                track_line = create_track_lines(sample_id, haplotype)
                outfile.write(f"{track_line}\n")
                
                # Read and process data
                df = pd.read_csv(
                    bed_file, 
                    sep='\t', 
                    header=None,
                    skiprows=1,
                    names=['chrom', 'start', 'end', 'name'],
                    dtype={
                        'chrom': str,
                        'start': int,
                        'end': int,
                        'name': str
                    }
                )
                
                # Filter for specific haplotype
                mask = df['name'].str.startswith(f"{haplotype}_")
                hap_regions = df[mask].copy()
                
                if not hap_regions.empty:
                    # Add required columns for IGV
                    hap_regions['sample_id'] = sample_id
                    hap_regions['score'] = 1000
                    hap_regions['strand'] = '.'
                    hap_regions['color'] = hap_regions['name'].apply(methylation_status_to_rgb)
                    
                    # Stagger overlapping regions within this sample
                    hap_regions = stagger_overlaps(hap_regions)
                    
                    # Select and order columns for IGV BED format
                    igv_bed = hap_regions[['chrom', 'start', 'end', 'sample_id', 'score', 
                                         'strand', 'start', 'end', 'color']]
                    
                    # Write sample data
                    igv_bed.to_csv(outfile, sep='\t', header=False, index=False)
                    logging.info(f"Processed {len(hap_regions)} {haplotype} regions for sample {sample_id}")
                
            except Exception as e:
                logging.error(f"Error processing {sample_id}: {e}")
                continue
                
    logging.info(f"Saved IGV BED file to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Create IGV visualization files for methylation data")
    parser.add_argument("--input_dir", required=True, help="Input directory containing BED files")
    parser.add_argument("--output_dir", required=True, help="Output directory for IGV files")
    parser.add_argument("--sample_list", help="Optional file containing sample IDs to process (in desired order)")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read sample list if provided
    sample_list = None
    if args.sample_list:
        sample_list = read_sample_list(args.sample_list)
    
    # Process each haplotype
    for haplotype in ['H1', 'H2']:
        aggregate_haplotype_regions(args.input_dir, args.output_dir, haplotype, sample_list)

if __name__ == "__main__":
    main()