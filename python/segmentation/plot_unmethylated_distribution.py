import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import logging

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def load_methylation_data(file_path):
    try:
        # First check if file is empty
        with open(file_path) as f:
            first_line = f.readline().strip()
            if not first_line:
                logging.warning(f"{file_path} is empty")
                return None
        
        # Read the file, skipping any comment or header lines
        df = pd.read_csv(file_path, sep='\t',
                        comment='#',  # Skip any comment lines
                        names=['chrom', 'start', 'end', 'label', 'score', 'strand'],
                        dtype={'chrom': str, 'start': str, 'end': str,  # Read as strings first
                              'label': str, 'score': str, 'strand': str},
                        low_memory=False)
        
        # Convert numeric columns after loading
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        
        # Drop any rows where conversion failed
        df = df.dropna(subset=['start', 'end'])
                        
        # Ensure required columns exist
        required_cols = ['chrom', 'start', 'end']
        if not all(col in df.columns for col in required_cols):
            logging.error(f"Missing required columns in {file_path}")
            return None
            
        # Keep only the columns we need
        df = df[required_cols]
        df['size'] = df['end'] - df['start']
        return df
    except pd.errors.EmptyDataError:
        logging.warning(f"{file_path} is empty")
        return None
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return None

def plot_distribution(input_dir, output_dir, sample_name):
    dfs = []
    input_dir = Path(input_dir)
    
    for hap in ['H1', 'H2']:
        for state in ['M', 'U']:
            file_path = input_dir / f"{hap}_{state}" / f"{sample_name}_{hap}_{state}.bed"
            if file_path.exists():
                df = load_methylation_data(file_path)
                if df is not None:
                    df['haplotype'] = hap
                    df['state'] = state
                    dfs.append(df)
                    logging.info(f"Loaded {file_path}")

    if not dfs:
        logging.error(f"No data found for sample {sample_name}")
        return

    combined_df = pd.concat(dfs, ignore_index=True)
    create_plots(combined_df, output_dir, sample_name)

def create_plots(df, output_dir, sample_name):
    plt.figure(figsize=(15, 6))
    
    plt.subplot(1, 2, 1)
    sns.histplot(data=df[df['state'] == 'M'],
                x='size',
                hue='haplotype',
                multiple="layer",
                palette={'H1': '#FF0000', 'H2': '#FF6666'},
                alpha=0.5,
                log_scale=True,
                bins=50)
    plt.title(f"Methylated Regions - {sample_name}")
    plt.xlabel("Region Size (bp)")
    plt.ylabel("Count")
    plt.grid(True, alpha=0.2)
    
    plt.subplot(1, 2, 2)
    sns.histplot(data=df[df['state'] == 'U'],
                x='size',
                hue='haplotype',
                multiple="layer",
                palette={'H1': '#0000FF', 'H2': '#6666FF'},
                alpha=0.5,
                log_scale=True,
                bins=50)
    plt.title(f"Unmethylated Regions - {sample_name}")
    plt.xlabel("Region Size (bp)")
    plt.grid(True, alpha=0.2)
    
    plt.suptitle(f"Methylation Region Size Distribution - {sample_name}", y=1.05)
    plt.tight_layout()
    
    output_path = Path(output_dir) / f"{sample_name}_methylation_distribution.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved plot to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Plot distribution of methylated/unmethylated regions.")
    parser.add_argument("--input-dir", required=True, help="Path to the input directory")
    parser.add_argument("--output-dir", required=True, help="Path to the output directory")
    parser.add_argument("--sample-name", required=True, help="Sample name to process")
    
    args = parser.parse_args()
    setup_logging()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    plot_distribution(args.input_dir, args.output_dir, args.sample_name)

if __name__ == "__main__":
    main()
