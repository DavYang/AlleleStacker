#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_bed_files(base_dir, condition):
    path = Path(base_dir) / condition
    dfs = []
    print(f"Loading BED files for {condition} from {path}")
    
    for bed_file in path.glob('*.bed'):
        df = pd.read_csv(bed_file, sep='\t')
        df['sample'] = bed_file.stem
        dfs.append(df)
        
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    else:
        print(f"No BED files found for {condition}")
        return None

def filter_chromosomes(df):
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    filtered_df = df[df['chrom'].isin(valid_chroms)]
    print(f"Filtered chromosomes: {filtered_df['chrom'].unique()}")
    
    return filtered_df

def plot_sample_distributions(df, output_dir, condition):
    if df is None:
        print(f"No data to plot for {condition}.")
        return
        
    output_dir = Path(output_dir)
    df = filter_chromosomes(df)
    samples = df['sample'].unique()
    
    # Create condition-level plots directory
    condition_dir = output_dir / condition
    condition_dir.mkdir(exist_ok=True, parents=True)
    print(f"Created or using {condition} directory at {condition_dir}")
    
    # Create samples subdirectory
    samples_dir = condition_dir / 'samples'
    samples_dir.mkdir(exist_ok=True, parents=True)
    print(f"Created or using samples directory at {samples_dir}")
    
    for sample in samples:
        sample_data = df[df['sample'] == sample]
        chromosomes = sorted(sample_data['chrom'].unique(), 
                           key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else float('inf'))
        
        fig, axes = plt.subplots(len(chromosomes), 1, figsize=(15, len(chromosomes)*0.5), sharex=True)
        fig.suptitle(f'Region Distribution - {sample}')
        
        for idx, chrom in enumerate(chromosomes):
            chrom_data = sample_data[sample_data['chrom'] == chrom]
            positions = (chrom_data['start'] + chrom_data['end']) / 2
            ax = axes[idx] if len(chromosomes) > 1 else axes
            ax.hist(positions, bins=100, alpha=0.7)
            ax.set_ylabel(chrom)
            ax.set_yticks([])
            
        plt.xlabel('Chromosome Position')
        plt.tight_layout()
        plot_path = samples_dir / f'{sample}_chromosome_density.png'
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"Saved plot for {sample} at {plot_path}")
        plt.close(fig)

    # Plot summary statistics
    fig_summary = plt.figure()
    sns.boxplot(x="chrom", y="start", data=df)
    plt.title("Summary Statistics")
    plt.tight_layout()
    summary_plot_path = condition_dir / 'summary_stats.png'
    plt.savefig(summary_plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved summary plot at {summary_plot_path}")
    
    # Save summary statistics to CSV
    summary_csv_path = condition_dir / 'summary_stats.csv'
    df_summary = df.groupby('chrom').agg({'start': ['mean', 'median', 'std']})
    df_summary.to_csv(summary_csv_path)
    print(f"Saved summary stats to {summary_csv_path}")

def main(base_dir, output_dir):
    conditions = ['H1_M', 'H1_U', 'H2_M', 'H2_U']
    
    for condition in conditions:
        print(f"Processing condition: {condition}")
        df = load_bed_files(base_dir, condition)
        
        if df is not None:
            plot_sample_distributions(df, output_dir, condition)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    args = parser.parse_args()
    
    print(f"Starting with base directory: {args.base_dir} and output directory: {args.output_dir}")
    main(args.base_dir, args.output_dir)
