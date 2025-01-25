import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import logging

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_methylation_data(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        required_cols = ['chrom', 'start', 'end']
        if not all(col in df.columns for col in required_cols):
            logging.error(f"Missing required columns in {file_path}")
            return None
        df['size'] = df['end'] - df['start']
        return df
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return None

def plot_distribution(input_dir, output_dir, sample_name):
    dfs = []
    input_dir = Path(input_dir)
    sample_output_dir = Path(output_dir) / sample_name
    sample_output_dir.mkdir(parents=True, exist_ok=True)
    
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
    create_size_position_plot(combined_df, sample_output_dir, sample_name)

def create_size_position_plot(df, output_dir, sample_name):
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    df = df[df['chrom'].isin(valid_chroms)].copy()
    
    for state, title, colors in [('M', 'Methylated', ('#FF0000', '#FF6666')), 
                                ('U', 'Unmethylated', ('#0000FF', '#6666FF'))]:
        plt.figure(figsize=(20, 8))
        state_df = df[df['state'] == state]
        
        for hap, color in zip(['H1', 'H2'], colors):
            hap_df = state_df[state_df['haplotype'] == hap]
            plt.vlines(x=hap_df['start'], 
                      ymin=0,
                      ymax=hap_df['size'],
                      color=color,
                      alpha=0.5,
                      label=hap)
        
        plt.yscale('log')
        plt.title(f"{title} Region Sizes by Position - {sample_name}")
        plt.xlabel("Genomic Position")
        plt.ylabel("Region Size (bp)")
        plt.grid(True, alpha=0.2)
        plt.legend()
        
        output_path = output_dir / f"{sample_name}_{state}_size_position.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved size-position plot to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Plot methylation regions by size and position.")
    parser.add_argument("--input-dir", required=True, help="Path to input directory")
    parser.add_argument("--output-dir", required=True, help="Path to output directory")
    parser.add_argument("--sample-name", required=True, help="Sample name to process")
    
    args = parser.parse_args()
    setup_logging()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    plot_distribution(args.input_dir, args.output_dir, args.sample_name)

if __name__ == "__main__":
    main()