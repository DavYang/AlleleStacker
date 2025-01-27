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

    # Filter out alternate haplotypes
    autosomes_and_sex_chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    df_filtered = df[df['chrom'].isin(autosomes_and_sex_chromosomes)]

    samples = df_filtered['sample'].unique()
    
    # Chromosome density plot
    chromosomes = sorted(df_filtered['chrom'].unique(), 
                         key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else float('inf'))
    
    for sample in samples:
        fig, axes = plt.subplots(len(chromosomes), 1, figsize=(15, len(chromosomes)*1.5), sharex=True)
        fig.suptitle(f'Region Distribution Across Chromosomes - {condition} - Sample: {sample}')
        
        for idx, chrom in enumerate(chromosomes):
            chrom_data = df_filtered[(df_filtered['chrom'] == chrom) & (df_filtered['sample'] == sample)]
            positions = (chrom_data['start'] + chrom_data['end']) / 2
            axes[idx].hist(positions, bins=100, alpha=0.7)
            axes[idx].set_ylabel(chrom)

        plt.xlabel('Chromosome Position')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(output_dir / f'{condition}_{sample}_chromosome_density.png', dpi=300, bbox_inches='tight')
        plt.close()

    # Size distribution
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_filtered, x='sample', y='size', showfliers=False)
    plt.xticks(rotation=45, ha='right')
    plt.yscale('log')
    plt.title(f'Region Size Distribution - {condition}')
    plt.tight_layout()
    plt.savefig(output_dir / f'{condition}_size_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Summary statistics
    stats = df_filtered.groupby('sample').agg({
        'size': ['count', 'mean', 'median', 'std'],
        'chrom': lambda x: len(x.unique())
    }).round(2)
    stats.to_csv(output_dir / f'{condition}_summary_stats.csv')

def main(base_dir, output_dir):
    conditions = ['H1_M', 'H1', 'H2']
    
    for condition in conditions:
        df = load_bed_files(base_dir, condition)
        plot_distributions(df, output_dir, condition)

if __name__ == "__main__":
    import sys
    base_dir = sys.argv[1]
    output_dir = sys.argv[2]
    main(base_dir, output_dir)
