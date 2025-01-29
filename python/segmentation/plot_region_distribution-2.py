#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import sys

# Centromere positions (hg38)
CENTROMERES = {
    'chr1': 122026460, 'chr2': 92188146, 'chr3': 90772459, 
    'chr4': 49708101, 'chr5': 46485901, 'chr6': 58830166,
    'chr7': 58054331, 'chr8': 43838887, 'chr9': 47367679,
    'chr10': 39254935, 'chr11': 51644205, 'chr12': 34856694,
    'chr13': 16000000, 'chr14': 16000000, 'chr15': 17000000,
    'chr16': 36311159, 'chr17': 22813680, 'chr18': 15460898,
    'chr19': 24498981, 'chr20': 26436233, 'chr21': 11288129,
    'chr22': 14300000, 'chrX': 58632012, 'chrY': 10104553
}

# Chromosome sizes (hg38)
CHROM_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
    'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

def load_bed_files(base_dir, condition):
    """Load and combine BED files for a given condition"""
    print(f"\nProcessing condition: {condition}")
    path = Path(base_dir) / condition
    dfs = []
    
    if not path.exists():
        print(f"Warning: Path {path} does not exist")
        return None
        
    bed_files = list(path.glob('*.bed'))
    print(f"Found {len(bed_files)} BED files")
    
    for bed_file in bed_files:
        print(f"Loading {bed_file.name}...")
        try:
            df = pd.read_csv(bed_file, sep='\t')
            df['sample'] = bed_file.stem
            dfs.append(df)
        except Exception as e:
            print(f"Error loading {bed_file}: {str(e)}")
            continue
    
    if not dfs:
        print("No valid BED files found")
        return None
        
    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Total regions loaded: {len(combined_df)}")
    return combined_df

def filter_chromosomes(df):
    """Filter for standard chromosomes only"""
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    filtered_df = df[df['chrom'].isin(valid_chroms)]
    print(f"Filtered from {len(df)} to {len(filtered_df)} regions")
    return filtered_df

def calculate_density(regions, chrom_size, bin_size=1000000):
    """Calculate region density considering region sizes"""
    print(f"Calculating density with {len(regions)} regions")
    bins = np.arange(0, chrom_size + bin_size, bin_size)
    density = np.zeros(len(bins)-1)
    
    for idx, region in enumerate(regions.itertuples()):
        if idx % 1000 == 0:  # Progress update every 1000 regions
            sys.stdout.write(f"\rProcessing region {idx+1}/{len(regions)}")
            sys.stdout.flush()
            
        # Find which bins this region overlaps
        start_bin = np.searchsorted(bins, region.start) - 1
        end_bin = np.searchsorted(bins, region.end) - 1
        
        for bin_idx in range(start_bin, end_bin + 1):
            if bin_idx >= 0 and bin_idx < len(density):
                # Calculate overlap with this bin
                bin_start = bins[bin_idx]
                bin_end = bins[bin_idx + 1]
                overlap_start = max(region.start, bin_start)
                overlap_end = min(region.end, bin_end)
                overlap = overlap_end - overlap_start
                density[bin_idx] += overlap / bin_size
    
    print("\nDensity calculation complete")
    return bins[:-1], density

def plot_sample_distributions(df, output_dir, condition):
    """Generate density plots for each sample"""
    if df is None:
        return
        
    print(f"\nGenerating plots for condition: {condition}")
    output_dir = Path(output_dir)
    df = filter_chromosomes(df)
    samples = df['sample'].unique()
    
    # Create output directories
    condition_dir = output_dir / condition
    condition_dir.mkdir(exist_ok=True)
    samples_dir = condition_dir / 'samples'
    samples_dir.mkdir(exist_ok=True)
    
    print(f"Processing {len(samples)} samples")
    
    for sample_idx, sample in enumerate(samples, 1):
        print(f"\nProcessing sample {sample_idx}/{len(samples)}: {sample}")
        sample_data = df[df['sample'] == sample]
        chromosomes = sorted(sample_data['chrom'].unique(), 
                           key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else float('inf'))
        
        fig, axes = plt.subplots(len(chromosomes), 1, figsize=(15, len(chromosomes)*0.5), sharex=True)
        fig.suptitle(f'Region Distribution - {sample}')
        
        print(f"Plotting {len(chromosomes)} chromosomes")
        for idx, chrom in enumerate(chromosomes):
            print(f"Processing chromosome {chrom}")
            chrom_data = sample_data[sample_data['chrom'] == chrom]
            ax = axes[idx] if len(chromosomes) > 1 else axes
            
            # Calculate and plot density
            if chrom in CHROM_SIZES:
                bins, density = calculate_density(chrom_data, CHROM_SIZES[chrom])
                ax.fill_between(bins, density, alpha=0.7)
                
                # Add centromere
                if chrom in CENTROMERES:
                    ax.axvline(x=CENTROMERES[chrom], color='red', linestyle='--', alpha=0.5, label='Centromere')
                
                # Add chromosome bands
                ax.set_xlim(0, CHROM_SIZES[chrom])
                
            ax.set_ylabel(chrom)
            ax.set_yticks([])
            ax.ticklabel_format(style='plain', axis='x')
            ax.grid(True, alpha=0.3)
            
            # Add coverage percentage
            total_covered = chrom_data['size'].sum()
            coverage_pct = (total_covered / CHROM_SIZES[chrom]) * 100 if chrom in CHROM_SIZES else 0
            ax.text(0.99, 0.8, f'{coverage_pct:.1f}% covered', 
                   transform=ax.transAxes, ha='right', fontsize=8)
        
        plt.xlabel('Chromosome Position (bp)')
        plt.tight_layout()
        
        output_file = samples_dir / f'{sample}_chromosome_density.png'
        print(f"Saving plot to {output_file}")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    print("\nGenerating summary statistics...")
    # Summary statistics
    summary_data = df.groupby(['chrom', 'sample'])['start'].count().reset_index()
    plt.figure(figsize=(12, 8))
    
    chromosomes = sorted(df['chrom'].unique(), 
                        key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else float('inf'))
    
    plt.boxplot([summary_data[summary_data['chrom'] == chrom]['start'] 
                for chrom in chromosomes])
    
    plt.xticks(range(1, len(chromosomes) + 1), chromosomes, rotation=90)
    plt.ylabel('Number of Regions')
    plt.title(f'Distribution of Region Counts by Chromosome - {condition}')
    plt.gca().yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    plt.tight_layout()
    output_file = condition_dir / 'summary_stats.png'
    print(f"Saving summary plot to {output_file}")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save summary statistics to CSV
    stats = df.groupby('sample').agg({
        'size': ['count', 'mean', 'median', 'std'],
        'chrom': lambda x: len(x.unique())
    }).round(2)
    output_file = condition_dir / 'summary_stats.csv'
    print(f"Saving summary statistics to {output_file}")
    stats.to_csv(output_file)

def main(base_dir, output_dir):
    print(f"Starting analysis...")
    print(f"Base directory: {base_dir}")
    print(f"Output directory: {output_dir}")
    
    conditions = ['H1_M', 'H1_U', 'H2_M', 'H2_U']
    print(f"Processing {len(conditions)} conditions: {', '.join(conditions)}")
    
    for condition in conditions:
        df = load_bed_files(base_dir, condition)
        plot_sample_distributions(df, output_dir, condition)
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', required=True, help='Base directory containing condition folders')
    parser.add_argument('--output_dir', required=True, help='Output directory for plots')
    args = parser.parse_args()
    main(args.base_dir, args.output_dir)
