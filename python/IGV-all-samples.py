import os
import sys
import glob
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
from tqdm import tqdm
import time
import subprocess

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def read_sample_list(sample_file):
    """Read and validate sample list from file"""
    try:
        with open(sample_file) as f:
            samples = [line.strip() for line in f if line.strip()]
        logging.info(f"Read {len(samples)} samples from {sample_file}")
        return samples
    except Exception as e:
        logging.error(f"Error reading sample list file {sample_file}: {e}")
        sys.exit(1)

def process_file(args):
    """Process a single file"""
    file_path, haplotype = args
    sample_id = Path(file_path).stem.replace('.meth_regions', '')
    
    try:
        # Read data
        df = pd.read_csv(
            file_path, 
            sep='\t', 
            header=None,
            skiprows=1,
            usecols=[0,1,2,3],
            names=['chrom', 'start_pos', 'end_pos', 'name'],
            dtype={
                'chrom': 'category',
                'start_pos': np.int32,
                'end_pos': np.int32,
                'name': str
            }
        )
        
        # Filter for haplotype
        mask = df['name'].str.startswith(f"{haplotype}_")
        if not mask.any():
            return None, sample_id, 0
            
        df = df[mask]
        region_count = len(df)
        
        # Add color based on methylation status
        df['color'] = np.where(df['name'].str.contains('_M'), '255,0,0', '0,0,255')
        df['sample_id'] = sample_id
        df['score'] = 1000
        df['strand'] = '.'
        
        return df, sample_id, region_count
        
    except Exception as e:
        logging.error(f"Error processing {file_path}: {e}")
        return None, sample_id, 0

def write_output_files(output_base, sample_data, ordered_samples, haplotype):
    """Write output files with consistent track positions"""
    igv_file = f"{output_base}/igv_cohort_{haplotype.lower()}.igv.bed"
    bed_file = f"{output_base}/igv_cohort_{haplotype.lower()}.bed"
    regions_written = 0
    
    # Write IGV file with tracks
    with open(igv_file, 'w') as igv_out, open(bed_file, 'w') as bed_out:
        # Write track definitions at the start of IGV file
        for i, sample_id in enumerate(ordered_samples):
            if sample_id in sample_data:
                track_line = (
                    f'track name="{sample_id}_{haplotype}" '
                    f'description="{sample_id} {haplotype} Methylation" '
                    f'visibility=full '
                    f'autoScale=off '
                    f'viewLimits=0:1000 '
                    f'useScore=0 '
                    f'itemRgb=On\n'
                )
                igv_out.write(track_line)

        # Process each chromosome
        all_chroms = set()
        for data in sample_data.values():
            all_chroms.update(data['chrom'].unique())

        for chrom in sorted(all_chroms):
            # For each chromosome, write samples in order
            for sample_id in ordered_samples:
                if sample_id in sample_data:
                    data = sample_data[sample_id]
                    chrom_data = data[data['chrom'] == chrom]
                    
                    if not chrom_data.empty:
                        # Sort by position within this sample/chromosome
                        chrom_data = chrom_data.sort_values(['start_pos', 'end_pos'])
                        
                        # Write regions to both files
                        for _, row in chrom_data.iterrows():
                            bed_line = [
                                row['chrom'],
                                str(row['start_pos']),
                                str(row['end_pos']),
                                f"{sample_id}_{haplotype}",
                                str(row['score']),
                                row['strand'],
                                str(row['start_pos']),
                                str(row['end_pos']),
                                row['color']
                            ]
                            bed_string = '\t'.join(bed_line) + '\n'
                            igv_out.write(bed_string)
                            bed_out.write(bed_string)
                            regions_written += 1
    
    # Create tabix index
    try:
        logging.info("Creating tabix index...")
        sort_cmd = f"sort -k1,1 -k2,2n {bed_file} > {bed_file}.sorted"
        subprocess.run(sort_cmd, shell=True, check=True)
        
        # Replace original with sorted file
        os.rename(f"{bed_file}.sorted", bed_file)
        
        # Compress with bgzip
        bgzip_cmd = f"bgzip -f {bed_file}"
        subprocess.run(bgzip_cmd, shell=True, check=True)
        
        # Index with tabix
        tabix_cmd = f"tabix -p bed {bed_file}.gz"
        subprocess.run(tabix_cmd, shell=True, check=True)
        
        logging.info("Successfully created tabix index")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during indexing: {e}")
    except Exception as e:
        logging.error(f"Unexpected error during indexing: {e}")
    
    return igv_file, f"{bed_file}.gz", regions_written

def process_single_haplotype():
    parser = argparse.ArgumentParser(description="Process single haplotype for IGV visualization")
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--haplotype", required=True, choices=['H1', 'H2'])
    parser.add_argument("--sample_list", required=True, help="File containing sample IDs in desired order")
    parser.add_argument("--threads", type=int, default=4)
    
    args = parser.parse_args()
    
    start_time = time.time()
    logging.info(f"Starting processing for {args.haplotype} haplotype")
    logging.info(f"Input directory: {args.input_dir}")
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Sample list file: {args.sample_list}")
    logging.info(f"Using {args.threads} threads")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read sample list
    ordered_samples = read_sample_list(args.sample_list)
    
    # Validate input files
    input_files = []
    missing_files = []
    for sample in ordered_samples:
        file_path = Path(args.input_dir) / f"{sample}.meth_regions.bed"
        if file_path.exists():
            input_files.append(file_path)
        else:
            missing_files.append(sample)
    
    if missing_files:
        logging.warning(f"Files not found for {len(missing_files)} samples: {', '.join(missing_files)}")
    
    if not input_files:
        logging.error("No valid input files found")
        sys.exit(1)
    
    logging.info(f"Found {len(input_files)} valid input files out of {len(ordered_samples)} samples")
    
    # Process files in parallel
    tasks = [(str(f), args.haplotype) for f in input_files]
    sample_data = {}
    total_regions = 0
    processed_samples = 0
    
    logging.info(f"Processing {len(tasks)} files in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_file, task): task for task in tasks}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing samples"):
            file_path, _ = futures[future]
            try:
                data, sample_id, region_count = future.result()
                if data is not None and not data.empty:
                    sample_data[sample_id] = data
                    total_regions += region_count
                    processed_samples += 1
                    logging.info(f"Processed {sample_id}: found {region_count} {args.haplotype} regions")
            except Exception as e:
                logging.error(f"Error processing {file_path}: {e}")
    
    if not sample_data:
        logging.error(f"No {args.haplotype} regions found in any samples")
        sys.exit(1)
    
    logging.info(f"Writing output files with {total_regions} total regions from {processed_samples} samples...")
    
    igv_file, indexed_file, regions_written = write_output_files(
        args.output_dir, sample_data, ordered_samples, args.haplotype
    )
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    # Final summary
    logging.info("=== Processing Summary ===")
    logging.info(f"Total processing time: {processing_time:.2f} seconds")
    logging.info(f"Samples processed: {processed_samples}/{len(ordered_samples)}")
    logging.info(f"Total regions found: {total_regions}")
    logging.info(f"Regions written: {regions_written}")
    logging.info(f"Average regions per sample: {total_regions/processed_samples:.1f}")
    logging.info(f"Processing speed: {total_regions/processing_time:.1f} regions/second")
    logging.info(f"IGV visualization file: {igv_file}")
    logging.info(f"Indexed BED file: {indexed_file}")
    logging.info("=== Complete ===")

if __name__ == "__main__":
    process_single_haplotype()
