import os
import pandas as pd
import sys
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def combine_meth_regions(input_dir, output_dir):
    """
    Combine methylation region files from H1_M, H1_U, H2_M, H2_U directories
    into IGV-compatible BED files
    """
    input_dir = Path(input_dir)
    logging.info(f"Processing input directory: {input_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define IGV-compatible RGB colors
    COLOR_MAP = {
        'M': '255,0,0',    # Red for methylated
        'U': '0,0,255'     # Blue for unmethylated
    }

    # Get all unique sample IDs
    sample_ids = set()
    for subdir in ['H1_M', 'H1_U', 'H2_M', 'H2_U']:
        dir_path = input_dir / subdir
        logging.info(f"Scanning directory: {dir_path}")
        
        if not dir_path.exists():
            logging.error(f"Directory not found: {dir_path}")
            continue
            
        for file_path in dir_path.glob('*.bed'):
            logging.info(f"Found file: {file_path.name}")
            # Extract sample ID - it's the first part before _H1_M, _H1_U, etc.
            sample_id = file_path.name.split('_')[0]
            sample_ids.add(sample_id)
            logging.info(f"Added sample ID: {sample_id}")

    logging.info(f"Found {len(sample_ids)} unique samples: {', '.join(sorted(sample_ids))}")

    # Process each sample
    for sample_id in sorted(sample_ids):
        for haplotype in ['H1', 'H2']:
            dfs = []
            
            # Read methylated regions
            m_file = input_dir / f"{haplotype}_M" / f"{sample_id}_{haplotype}_M.bed"
            if m_file.exists():
                logging.info(f"Reading methylated file: {m_file}")
                try:
                    m_df = pd.read_csv(m_file, sep='\t')
                    m_df['rgb'] = COLOR_MAP['M']
                    dfs.append(m_df)
                    logging.info(f"Read {len(m_df)} methylated regions")
                except Exception as e:
                    logging.error(f"Error reading {m_file}: {e}")

            # Read unmethylated regions
            u_file = input_dir / f"{haplotype}_U" / f"{sample_id}_{haplotype}_U.bed"
            if u_file.exists():
                logging.info(f"Reading unmethylated file: {u_file}")
                try:
                    u_df = pd.read_csv(u_file, sep='\t')
                    u_df['rgb'] = COLOR_MAP['U']
                    dfs.append(u_df)
                    logging.info(f"Read {len(u_df)} unmethylated regions")
                except Exception as e:
                    logging.error(f"Error reading {u_file}: {e}")

            if dfs:
                # Combine and sort regions
                combined_df = pd.concat(dfs, ignore_index=True)
                combined_df = combined_df.sort_values(['chrom', 'start'])

                # Create IGV-compatible BED format
                output_df = pd.DataFrame({
                    'chrom': combined_df['chrom'],
                    'start': combined_df['start'],
                    'end': combined_df['end'],
                    'name': f"{sample_id}_{haplotype}",
                    'score': 1000,
                    'strand': '.',
                    'thickStart': combined_df['start'],
                    'thickEnd': combined_df['end'],
                    'rgb': combined_df['rgb']
                })

                # Write output file
                output_file = Path(output_dir) / f"{sample_id}_{haplotype}.bed"
                output_df.to_csv(output_file, sep='\t', index=False, header=False)
                logging.info(f"Processed {sample_id} {haplotype}: {len(output_df)} regions")

            else:
                logging.warning(f"No data found for {sample_id} {haplotype}")

    logging.info("Processing complete. Check the output directory for results.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_meth_regions.py <input_dir> <output_dir>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    combine_meth_regions(input_directory, output_directory)