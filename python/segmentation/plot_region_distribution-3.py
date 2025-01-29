#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import matplotlib.gridspec as gridspec

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
            sample_name = bed_file.stem.split('_')[0]
            df['sample'] = sample_name
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

def format_genomic_position(x, p):
    """Format genomic positions in Mb"""
    return f'{int(x/1e6)}'

def plot_chromosome_distribution(fig, position, chrom, data, is_methylated=True, is_last=False):
    """Create a histogram plot for each chromosome"""
    ax = fig.add_subplot(position)
    
    # Plot data
    if len(data) > 0:
        positions = (data['start'] + data['end']) / 2
        color = '#FF6B6B' if is_methylated else '#4ECDC4'
        ax.hist(positions, bins=500, alpha=0.6, color=color,
               range=(0, CHROM_SIZES[chrom]))
    
    # Remove y-axis elements
    ax.set_yticks([])
    ax.set_yticklabels([])
    
    # Create a clear separation between plots
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Handle x-axis
    if not is_last:
        ax.set_xticklabels([])
        ax.spines['bottom'].set_visible(False)
    else:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_genomic_position))
        ax.set_xlabel('Position (Mb)')
    
    # Set x limits based on individual chromosome size
    ax.set_xlim(-0.01 * CHROM_SIZES[chrom], 1.01 * CHROM_SIZES[chrom])
    
    # Add centromere line
    if chrom in CENTROMERES:
        ax.axvline(x=CENTROMERES[chrom], color='black', linestyle='--',
                  linewidth=1, alpha=0.7)
    
    # Add subtle grid
    ax.grid(True, alpha=0.1, linestyle=':')
    
    # Make chromosome label prominent
    ax.text(-0.12, 0.5, chrom, transform=ax.transAxes,
            fontsize=10, fontweight='bold', ha='right', va='center')
    
    # Add statistics in a more subtle way
    count = len(data)
    size = int(data['size'].mean()) if len(data) > 0 else 0
    stats = f'{count:,} regions, {size:,}bp mean'
    ax.text(1.02, 0.5, stats, transform=ax.transAxes,
            fontsize=8, alpha=0.7, va='center')

def plot_methylation_comparison(base_dir, output_dir, haplotype):
    """Generate separate methylation and unmethylation plots"""
    
    # Set style
    plt.style.use('default')
    sns.set_theme(style="white")
    
    def load_and_process_data(condition):
        df = load_bed_files(base_dir, f"{haplotype}_{condition}")
        if df is not None:
            df = filter_chromosomes(df)
        return df
    
    # Load methylated and unmethylated data
    m_data = load_and_process_data('M')
    u_data = load_and_process_data('U')
    
    if m_data is None or u_data is None:
        print(f"Error: Could not load data for haplotype {haplotype}")
        return
    
    # Get unique samples
    samples = sorted(set([name.split('_')[0] for name in 
                        pd.concat([pd.Series(m_data['sample'].unique()), 
                                 pd.Series(u_data['sample'].unique())]).unique()]))
    
    for sample_idx, sample in enumerate(samples, 1):
        print(f"\nProcessing sample {sample_idx}/{len(samples)}: {sample}")
        
        # Get data for this sample
        m_sample = m_data[m_data['sample'] == sample]
        u_sample = u_data[u_data['sample'] == sample]
        
        # Get chromosomes
        chrom_list = list(set(list(m_sample['chrom'].unique()) + 
                            list(u_sample['chrom'].unique())))
        chromosomes = sorted(chrom_list, 
                           key=lambda x: int(x.replace('chr', '')) 
                           if x.replace('chr', '').isdigit() else float('inf'))
        
        # Create separate figures for methylated and unmethylated data
        for is_methylated in [True, False]:
            data = m_sample if is_methylated else u_sample
            title_prefix = 'Methylated' if is_methylated else 'Unmethylated'
            
            # Create figure with more horizontal space for labels
            fig = plt.figure(figsize=(15, len(chromosomes) * 0.6))
            gs = gridspec.GridSpec(len(chromosomes), 1)
            gs.update(hspace=0.1)  # Small space between plots for better readability
            
            # Add title
            fig.suptitle(f'{title_prefix} Region Distribution - {sample} ({haplotype})', 
                        y=0.95, fontsize=14, fontweight='bold')
            
            # Process each chromosome
            for idx, chrom in enumerate(chromosomes):
                chrom_data = data[data['chrom'] == chrom]
                is_last = idx == len(chromosomes) - 1
                plot_chromosome_distribution(fig, gs[idx], chrom, chrom_data, 
                                          is_methylated, is_last)
            
            # Adjust layout to ensure labels are visible
            plt.subplots_adjust(left=0.15, right=0.85)
            
            # Save plot
            output_dir = Path(output_dir)
            haplotype_dir = output_dir / haplotype
            haplotype_dir.mkdir(parents=True, exist_ok=True)
            
            prefix = 'methylated' if is_methylated else 'unmethylated'
            output_file = haplotype_dir / f'{sample}_{prefix}_distribution.png'
            print(f"Saving plot to {output_file}")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate methylation distribution plots')
    parser.add_argument('--base_dir', required=True, help='Base directory containing condition folders')
    parser.add_argument('--output_dir', required=True, help='Output directory for plots')
    args = parser.parse_args()
    
    print("Starting methylation visualization...")
    print(f"Base directory: {args.base_dir}")
    print(f"Output directory: {args.output_dir}")
    
    for haplotype in ['H1', 'H2']:
        plot_methylation_comparison(args.base_dir, args.output_dir, haplotype)
    
    print("\nVisualization complete!")

if __name__ == "__main__":
    main()