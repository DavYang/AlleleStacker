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
        df = pd.read_csv(file_path, sep='\t')
        required_cols = ['chrom', 'start', 'end']
        if not all(col in df.columns for col in required_cols):
            logging.error(f"Missing required columns in {file_path}")
            return None
            
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
    # Calculate cumulative chromosome lengths for linearization
    CHROMOSOME_LENGTHS = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
        'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
        'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
        'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
        'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
        'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }
    
    # Calculate cumulative positions
    cumulative_pos = {}
    current_pos = 0
    for chrom, length in chrom_lengths.items():  # Use chrom_lengths here
        cumulative_pos[chrom] = current_pos
        current_pos += length

    # Function to linearize position
    def get_linear_pos(row):
        return cumulative_pos[row['chrom']] + row['start']

    # Add linearized position
    df['linear_pos'] = df.apply(get_linear_pos, axis=1)
    
    # Create separate plots for methylated and unmethylated
    for state, title, colors in [('M', 'Methylated', ('#FF0000', '#FF6666')), 
                                ('U', 'Unmethylated', ('#0000FF', '#6666FF'))]:
        plt.figure(figsize=(20, 8))
        state_df = df[df['state'] == state]
        
        # Plot regions as lines to show their size
        for hap, color in zip(['H1', 'H2'], colors):
            hap_df = state_df[state_df['haplotype'] == hap]
            for _, row in hap_df.iterrows():
                plt.hlines(y=row['size'], 
                         xmin=row['linear_pos'],
                         xmax=row['linear_pos'] + row['size'],
                         color=color,
                         alpha=0.5,
                         linewidth=1)
        
        # Add chromosome boundary lines and labels
        for chrom, pos in cumulative_pos.items():
            plt.axvline(x=pos, color='gray', linestyle='--', alpha=0.3)
            plt.text(pos + chrom_lengths[chrom] / 2, # Use chrom_lengths here
                    plt.ylim()[1], 
                    chrom.replace('chr', ''),
                    ha='center',
                    va='bottom')
        
        plt.title(f"{title} Regions by Chromosome - {sample_name}")
        plt.xlabel("Genomic Position")
        plt.ylabel("Region Size (bp)")
        plt.yscale('log')
        plt.grid(True, alpha=0.2)
        plt.legend(['H1', 'H2'])
        
        output_path = output_dir / f"{sample_name}_{state}_linear_chromosome_dist.png"
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
