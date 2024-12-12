import os
import sys
import glob
import logging
import pandas as pd
import pybedtools

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def read_sample_list(sample_file):
    """Read sample IDs from a text file"""
    try:
        with open(sample_file, 'r') as f:
            # Read lines and remove whitespace
            samples = [line.strip() for line in f if line.strip()]
        logging.info(f"Loaded {len(samples)} samples from {sample_file}")
        return set(samples)  # Convert to set for faster lookup
    except Exception as e:
        logging.error(f"Error reading sample file {sample_file}: {e}")
        sys.exit(1)

def get_sample_id(filename):
    """Extract sample ID from the filename"""
    basename = os.path.basename(filename)
    sample_id = basename.replace('.meth_regions.bed', '')
    return sample_id

def aggregate_haplotype_regions(input_dir: str, output_dir: str, haplotype: str, sample_list: set = None):
    """
    Aggregate all regions for a specific haplotype (H1 or H2) from all samples into a single BED file.
    Includes sample ID in the output.
    
    Args:
        input_dir: Directory containing *meth_regions.bed files
        output_dir: Directory where output files will be saved
        haplotype: Either 'H1' or 'H2'
        sample_list: Optional set of sample IDs to process
    """
    logging.info(f"Processing {haplotype} regions...")
    
    # Output file path
    output_file = os.path.join(output_dir, f"{haplotype}_regions.bed")
    
    # Get all meth_regions.bed files
    bed_files = glob.glob(os.path.join(input_dir, "*meth_regions.bed"))
    if not bed_files:
        logging.error(f"No *meth_regions.bed files found in {input_dir}")
        return
    
    # Read and combine all haplotype-specific regions
    all_regions = []
    processed_samples = 0
    skipped_samples = 0
    
    for bed_file in bed_files:
        # Get sample ID from filename
        sample_id = get_sample_id(bed_file)
        
        # Skip if not in sample list
        if sample_list is not None and sample_id not in sample_list:
            skipped_samples += 1
            continue
        
        try:
            # Read BED file with specified dtypes, skipping header
            df = pd.read_csv(
                bed_file, 
                sep='\t', 
                header=None,
                skiprows=1,  # Skip the header row
                names=['chrom', 'start', 'end', 'name'],
                dtype={
                    'chrom': str,
                    'start': int,
                    'end': int,
                    'name': str
                }
            )
            
            # Filter for specific haplotype and create new DataFrame
            mask = df['name'].str.startswith(f"{haplotype}_")
            hap_regions = df[mask].copy()
            
            if not hap_regions.empty:
                # Add sample ID column using loc
                hap_regions.loc[:, 'sample_id'] = sample_id
                all_regions.append(hap_regions)
                processed_samples += 1
                logging.info(f"Found {len(hap_regions)} {haplotype} regions in sample {os.path.basename(bed_file)}")
        
        except Exception as e:
            logging.error(f"Error processing {bed_file}: {e}")
            continue
    
    # Log processing summary
    logging.info(f"Processed {processed_samples} samples")
    if sample_list is not None:
        logging.info(f"Skipped {skipped_samples} samples not in sample list")
    
    if not all_regions:
        logging.error(f"No {haplotype} regions found in any files")
        return
    
    # Combine all regions
    combined_regions = pd.concat(all_regions, ignore_index=True)
    logging.info(f"Total {haplotype} regions: {len(combined_regions)}")
    
    # Select and order columns
    combined_regions = combined_regions[['chrom', 'start', 'end', 'name', 'sample_id']]
    
    # Convert to BedTool, sort, and save
    try:
        bed = pybedtools.BedTool.from_dataframe(combined_regions)
        sorted_bed = bed.sort()
        
        # Save sorted BED file
        sorted_bed.saveas(output_file)
        logging.info(f"Saved sorted regions to {output_file}")
        
    except Exception as e:
        logging.error(f"Error during BED processing: {e}")
        return

def main():
    if len(sys.argv) not in [3, 4]:
        print("Usage: python aggregate_haps.py <input_dir> <output_dir> [sample_list.txt]")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Read sample list if provided
    sample_list = None
    if len(sys.argv) == 4:
        sample_file = sys.argv[3]
        sample_list = read_sample_list(sample_file)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each haplotype
    for haplotype in ['H1', 'H2']:
        aggregate_haplotype_regions(input_dir, output_dir, haplotype, sample_list)

if __name__ == "__main__":
    main()