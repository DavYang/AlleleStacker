#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

def load_bed_files(base_dir, condition):
    path = Path(base_dir) / condition
    dfs = []
    
    for bed_file in path.glob('*.bed'):
        df = pd.read_csv(bed_file, sep='\t')
        df['sample'] = bed_file.stem
        dfs.append(df)
    
    return pd.concat(dfs, ignore_index=True) if dfs else None

def plot_distributions(df, output_dir, condition):
    if df is None:
        return
        
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Chromosome density plot
    chromosomes = sorted(df['chrom'].unique(), 
                        key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else float('inf'))
    
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(15, len(chromosomes)*1.5), sharex=True)
    fig.suptitle(f'Region Distribution Across Chromosomes - {condition}')
    
    for idx, chrom in enumerate(chromosomes):
        chrom_data = df[df['chrom'] == chrom]
        positions = (chrom_data['start'] + chrom_data['end']) / 2
        axes[idx].hist(positions, bins=100, alpha=0.7)
        axes[idx].set_ylabel(chrom)
        
    plt.xlabel('Chromosome Position')
    plt.tight_layout()
    plt.savefig(output_dir / f'{condition}_chromosome_density.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Size distribution
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x='sample', y='size', showfliers=False)
    plt.xticks(rotation=45, ha='right')
    plt.yscale('log')
    plt.title(f'Region Size Distribution - {condition}')
    plt.tight_layout()
    plt.savefig(output_dir / f'{condition}_size_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Summary statistics
    stats = df.groupby('sample').agg({
        'size': ['count', 'mean', 'median', 'std'],
        'chrom': lambda x: len(x.unique())
    }).round(2)
    stats.to_csv(output_dir / f'{condition}_summary_stats.csv')

def main(base_dir, output_dir):
    conditions = ['H1_M', 'H1_U', 'H2_M', 'H2_U']
    
    for condition in conditions:
        df = load_bed_files(base_dir, condition)
        plot_distributions(df, output_dir, condition)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    
    args = parser.parse_args()
    main(args.base_dir, args.output_dir)