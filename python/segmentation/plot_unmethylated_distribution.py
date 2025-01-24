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
        # First check file content
        with open(file_path) as f:
            lines = f.readlines()
            if not lines:
                logging.warning(f"{file_path} is empty")
                return None
            
            # Skip header lines until we find valid data
            data_lines = []
            for line in lines:
                if line.strip() and not line.startswith('#'):
                    fields = line.strip().split('\t')
                    try:
                        # Test if second and third fields can be converted to int
                        int(fields[1])
                        int(fields[2])
                        data_lines.append(line)
                    except (ValueError, IndexError):
                        continue
                        
        if not data_lines:
            logging.warning(f"No valid data found in {file_path}")
            return None
            
        # Create DataFrame from valid lines
        df = pd.read_csv(pd.io.common.StringIO(''.join(data_lines)), 
                        sep='\t',
                        names=['chrom', 'start', 'end', 'label', 'score', 'strand'],
                        dtype={'chrom': str, 'start': int, 'end': int,
                              'label': str, 'score': str, 'strand': str})
        
        # Convert numeric columns after loading
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        
        # Drop any rows where conversion failed
        # Ensure required columns exist
        required_cols = ['chrom', 'start', 'end']
        if not all(col in df.columns for col in required_cols):
            logging.error(f"Missing required columns in {file_path}")
            return None

        # Keep only the columns we need
        df = df[required_cols]

        # Convert to numeric and drop rows with NaN in 'start' and 'end'
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        df = df.dropna(subset=['start', 'end'])
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
    # Create sample-specific output directory
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
    create_plots(combined_df, sample_output_dir, sample_name)
    create_chromosome_plot(combined_df, sample_output_dir, sample_name)


def create_chromosome_plot(df, output_dir, sample_name):
    plt.figure(figsize=(15, 8))
    
    # Define chromosome order
    chrom_order = ([str(i) for i in range(1, 23)] + ['X', 'Y'])
    df['chrom'] = pd.Categorical(df['chrom'].str.replace('chr', ''), 
                                categories=chrom_order, 
                                ordered=True)
    
    # Create separate plots for methylated and unmethylated
    for state, title, colors in [('M', 'Methylated', ('#FF0000', '#FF6666')), 
                                ('U', 'Unmethylated', ('#0000FF', '#6666FF'))]:
        plt.figure(figsize=(15, 6))
        state_df = df[df['state'] == state]
        
        for hap, color in zip(['H1', 'H2'], colors):
            hap_df = state_df[state_df['haplotype'] == hap]
            plt.scatter(hap_df['chrom'], 
                       hap_df['start'], 
                       c=color, 
                       alpha=0.5, 
                       s=10, 
                       label=f'{hap}')
        
        plt.title(f"{title} Regions by Chromosome - {sample_name}")
        plt.xlabel("Chromosome")
        plt.ylabel("Position (bp)")
        plt.yscale('log')
        plt.grid(True, alpha=0.2)
        plt.legend()
        plt.xticks(rotation=45)
        
        output_path = output_dir / f"{sample_name}_{state}_chromosome_dist.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved chromosome plot to {output_path}")

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
    
    output_path = output_dir / f"{sample_name}_methylation_distribution.png"
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
