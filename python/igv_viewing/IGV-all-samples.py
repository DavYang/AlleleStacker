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
    """Process methylated and unmethylated files for a sample"""
    sample_id, input_dir, haplotype = args
    
    try:
        # Read methylated and unmethylated files
        m_file = Path(input_dir) / f"H{haplotype}_M" / f"{sample_id}_H{haplotype}_M.bed"
        u_file = Path(input_dir) / f"H{haplotype}_U" / f"{sample_id}_H{haplotype}_U.bed"
        
        dfs = []
        for file_path, is_meth in [(m_file, True), (u_file, False)]:
            if file_path.exists():
                df = pd.read_csv(
                    file_path,
                    sep='\t',
                    dtype={
                        'chrom': 'category',
                        'start': np.int32,
                        'end': np.int32,
                        'summary_label': str,
                        'size': np.int32
                    }
                )
                df['color'] = '255,0,0' if is_meth else '0,0,255'
                dfs.append(df)
        
        if not dfs:
            return None, sample_id, 0
            
        # Combine methylated and unmethylated data
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df = combined_df.sort_values(['chrom', 'start'])
        
        # Add required columns for IGV
        combined_df['sample_id'] = sample_id
        combined_df['score'] = 1000
        combined_df['strand'] = '.'
        
        return combined_df, sample_id, len(combined_df)
        
    except Exception as e:
        logging.error(f"Error processing {sample_id}: {e}")
        return None, sample_id, 0

def write_output_files(output_base, sample_data, ordered_samples, haplotype):
    """Write output files with consistent track positions"""
    igv_file = f"{output_base}/igv_cohort_h{haplotype.lower()}.igv.bed"
    bed_file = f"{output_base}/igv_cohort_h{haplotype.lower()}.bed"
    regions_written = 0
    
    with open(igv_file, 'w') as igv_out, open(bed_file, 'w') as bed_out:
        # Write track definitions
        for sample_id in ordered_samples:
            if sample_id in sample_data:
                track_line = (
                    f'track name="{sample_id}_H{haplotype}" '
                    f'description="{sample_id} H{haplotype} Methylation" '
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
            for sample_id in ordered_samples:
                if sample_id in sample_data:
                    data = sample_data[sample_id]
                    chrom_data = data[data['chrom'] == chrom]
                    
                    if not chrom_data.empty:
                        for _, row in chrom_data.iterrows():
                            bed_line = [
                                row['chrom'],
                                str(row['start']),
                                str(row['end']),
                                f"{sample_id}_H{haplotype}",
                                str(row['score']),
                                row['strand'],
                                str(row['start']),
                                str(row['end']),
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
        os.rename(f"{bed_file}.sorted", bed_file)
        subprocess.run(f"bgzip -f {bed_file}", shell=True, check=True)
        subprocess.run(f"tabix -p bed {bed_file}.gz", shell=True, check=True)
        logging.info("Successfully created tabix index")
    except Exception as e:
        logging.error(f"Error during indexing: {e}")
    
    return igv_file, f"{bed_file}.gz", regions_written

def process_single_haplotype():
    parser = argparse.ArgumentParser(description="Process single haplotype for IGV visualization")
    parser.add_argument("--input_dir", required=True, help="Directory containing H1_M, H1_U, H2_M, H2_U subdirectories")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--haplotype", required=True, choices=['1', '2'])
    parser.add_argument("--sample_list", required=True, help="File containing sample IDs in desired order")
    parser.add_argument("--threads", type=int, default=4)
    
    args = parser.parse_args()
    
    start_time = time.time()
    logging.info(f"Starting processing for H{args.haplotype} haplotype")
    logging.info(f"Input directory: {args.input_dir}")
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Using {args.threads} threads")
    
    os.makedirs(args.output_dir, exist_ok=True)
    ordered_samples = read_sample_list(args.sample_list)
    
    # Validate input files
    missing_files = []
    for sample in ordered_samples:
        m_file = Path(args.input_dir) / f"H{args.haplotype}_M" / f"{sample}_H{args.haplotype}_M.bed"
        u_file = Path(args.input_dir) / f"H{args.haplotype}_U" / f"{sample}_H{args.haplotype}_U.bed"
        if not (m_file.exists() or u_file.exists()):
            missing_files.append(sample)
    
    if missing_files:
        logging.warning(f"No files found for {len(missing_files)} samples: {', '.join(missing_files)}")
        ordered_samples = [s for s in ordered_samples if s not in missing_files]
    
    if not ordered_samples:
        logging.error("No valid input files found")
        sys.exit(1)
    
    # Process files in parallel
    tasks = [(s, args.input_dir, args.haplotype) for s in ordered_samples]
    sample_data = {}
    total_regions = 0
    processed_samples = 0
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_file, task): task for task in tasks}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing samples"):
            sample_id, _, _ = futures[future]
            try:
                data, sample_id, region_count = future.result()
                if data is not None and not data.empty:
                    sample_data[sample_id] = data
                    total_regions += region_count
                    processed_samples += 1
                    logging.info(f"Processed {sample_id}: found {region_count} H{args.haplotype} regions")
            except Exception as e:
                logging.error(f"Error processing {sample_id}: {e}")
    
    if not sample_data:
        logging.error(f"No H{args.haplotype} regions found in any samples")
        sys.exit(1)
    
    logging.info(f"Writing output files with {total_regions} total regions from {processed_samples} samples...")
    
    igv_file, indexed_file, regions_written = write_output_files(
        args.output_dir, sample_data, ordered_samples, args.haplotype
    )
    
    end_time = time.time()
    processing_time = end_time - start_time
    
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