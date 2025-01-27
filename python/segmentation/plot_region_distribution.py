#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_bed_files(base_dir, condition):
    path = Path(base_dir) / condition
    dfs = []
    for bed_file in path.glob('*.bed'):
        df = pd.read_csv(bed_file, sep='\t')
        df['sample'] = bed_file.stem
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else None

def filter_chromosomes(df):
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    return df[df['chrom'].isin(valid_chroms)]

def plot_sample_distributions(df, output_dir, condition):
    if df is None:
        return
        
    output_dir = Path(output_dir)
    df = filter_chromosomes(df)
    samples = df['sample'].unique()
    
    # Create condition-level plots directory
    condition_dir = output_dir / condition
    condition_dir.mkdir(exist_ok=True)
    
    # Create samples subdirectory
    samples_dir = condition_dir / 'samples'
    samples_dir.mkdir(exist_ok=True)
    
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
        plt.savefig(samples_dir / f'{sample}_chromosome_density.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Condition-level plots
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x='sample', y='size', showfliers=False)
    plt.xticks(rotation=45, ha='right')
    plt.yscale('log')
    plt.title(f'Region Size Distribution - {condition}')
    plt.tight_layout()
    plt.savefig(condition_dir / 'size_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    stats = df.groupby('sample').agg({
        'size': ['count', 'mean', 'median', 'std'],
        'chrom': lambda x: len(x.unique())
    }).round(2)
    stats.to_csv(condition_dir / 'summary_stats.csv')

def main(base_dir, output_dir):
    conditions = ['H1_M', 'H1_U', 'H2_M', 'H2_U']
    for condition in conditions:
        df = load_bed_files(base_dir, condition)
        plot_sample_distributions(df, output_dir, condition)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    args = parser.parse_args()
    main(args.base_dir, args.output_dir)