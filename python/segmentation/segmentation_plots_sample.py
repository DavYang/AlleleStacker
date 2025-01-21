#!/usr/bin/env python3
"""
Methylation Region Analysis Script

Analyzes methylation regions from sample data and generates visualizations including:
- Region count plots
- Size distribution plots
- Chromosome distribution plots
- Summary statistics
"""

import logging
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')

# Constants
PLOT_STYLE = {
    'colors': {
        'Hap1-M': '#FF0000',  # Red
        'Hap1-U': '#0000FF',  # Blue
        'Hap2-M': '#FF6666',  # Light red
        'Hap2-U': '#6666FF'   # Light blue
    },
    'fontsize': {
        'title': 20,
        'label': 15,
        'tick': 13
    }
}

# Category definitions
CATEGORIES = ['Hap1-M', 'Hap1-U', 'Hap2-M', 'Hap2-U']

# Define chromosome order
CHROM_ORDER = ([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])

def read_methylation_file(file_path, haplotype):
    """
    Read and process a methylation data file
    
    Args:
        file_path: Path to the methylation BED file
        haplotype: Haplotype identifier (e.g., 'Hap1' or 'Hap2')
    
    Returns:
        DataFrame with processed methylation data or None if error
    """
    try:
        df = pd.read_csv(file_path, 
                        sep='\t',
                        comment='#',
                        names=['chrom', 'start', 'end', 'summary_label'],
                        usecols=['chrom', 'start', 'end', 'summary_label'])
        
        df['size'] = df['end'] - df['start']
        df['category'] = f"{haplotype}-" + df['summary_label']
        return df
        
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return None

def create_distribution_plot(data, category, sample_name, output_path):
    """
    Create size distribution plot for methylation regions
    
    Args:
        data: DataFrame containing size and category information
        category: Type of region (Methylated or Unmethylated)
        sample_name: Name of the sample
        output_path: Path to save the plot
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Calculate bin edges in log space
    bins = np.logspace(np.log10(data['size'].min()), 
                      np.log10(data['size'].max()), 
                      50)
    
    for i, cat in enumerate(data['category'].unique()):
        subset = data[data['category'] == cat]
        # Create step histogram with filled area
        counts, edges = np.histogram(subset['size'], bins=bins)
        plt.stairs(counts, edges, 
                  label=cat,
                  color=PLOT_STYLE['colors'][cat],
                  alpha=0.3,  # Fill transparency
                  fill=True)
        # Add line on top for better visibility
        plt.stairs(counts, edges,
                  color=PLOT_STYLE['colors'][cat],
                  linewidth=2,
                  alpha=1.0,
                  linestyle=['-', '--'][i])  # Different line styles for each haplotype
    
    plt.title(f'{category} Region Size Distribution - {sample_name}', 
              fontsize=PLOT_STYLE['fontsize']['title'])
    plt.xlabel('Region Size (bp)', fontsize=PLOT_STYLE['fontsize']['label'])
    plt.ylabel('Count', fontsize=PLOT_STYLE['fontsize']['label'])
    plt.xticks(fontsize=PLOT_STYLE['fontsize']['tick'])
    plt.yticks(fontsize=PLOT_STYLE['fontsize']['tick'])
    
    # Set log scales
    plt.xscale('log')
    plt.yscale('log')
    
    # Improve grid
    plt.grid(True, which="major", ls="-", alpha=0.2)
    plt.grid(True, which="minor", ls=":", alpha=0.1)
    
    # Improve legend
    plt.legend(frameon=True, framealpha=1.0, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_chromosome_distribution_plot(df, sample_name, output_path):
    """
    Create distribution plot of regions across chromosomes
    
    Args:
        df: DataFrame containing chromosome and category information
        sample_name: Name of the sample
        output_path: Path to save the plot
    """
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Get counts per chromosome for each category
    bar_width = 0.2
    for i, category in enumerate(CATEGORIES):
        # Get chromosome counts for this category
        subset = df[df['category'] == category]
        counts = subset['chrom'].value_counts().reindex(CHROM_ORDER, fill_value=0)
        
        # Plot bars with offset
        x = np.arange(len(CHROM_ORDER)) + (i * bar_width)
        plt.bar(x, counts.values, 
               width=bar_width,
               label=category,
               color=PLOT_STYLE['colors'][category],
               alpha=0.6)
    
    # Customize plot
    plt.title(f'Chromosome Distribution of Methylation Regions - {sample_name}', 
             fontsize=PLOT_STYLE['fontsize']['title'])
    plt.xlabel('Chromosome', fontsize=PLOT_STYLE['fontsize']['label'])
    plt.ylabel('Number of Regions', fontsize=PLOT_STYLE['fontsize']['label'])
    
    # Set x-axis ticks at center of each group
    group_center = bar_width * (len(CATEGORIES) - 1) / 2
    plt.xticks(np.arange(len(CHROM_ORDER)) + group_center, 
               [chrom.replace('chr', '') for chrom in CHROM_ORDER],
               rotation=45,
               ha='right',
               fontsize=PLOT_STYLE['fontsize']['tick'])
    
    # Add grid
    plt.grid(True, axis='y', linestyle='--', alpha=0.3)
    
    # Customize legend
    plt.legend(frameon=True, 
              framealpha=1.0, 
              loc='upper right',
              fontsize=PLOT_STYLE['fontsize']['tick'])
    
    plt.yscale('log')  # Use log scale for y-axis
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_sample_plots(df, sample_name, output_dir):
    """
    Create all sample-specific plots
    
    Args:
        df: DataFrame containing all sample data
        sample_name: Name of the sample
        output_dir: Directory to save plots
    """
    # Count plot
    fig, ax = plt.subplots(figsize=(10, 6))
    counts = df['category'].value_counts().reindex(CATEGORIES, fill_value=0)
    
    bars = ax.bar(range(len(CATEGORIES)), counts.values, 
                 color=[PLOT_STYLE['colors'][cat] for cat in CATEGORIES])
    
    ax.set_title(f'Methylation Region Counts - {sample_name}', 
                 fontsize=PLOT_STYLE['fontsize']['title'])
    ax.set_xlabel('Category', fontsize=PLOT_STYLE['fontsize']['label'])
    ax.set_ylabel('Count', fontsize=PLOT_STYLE['fontsize']['label'])
    ax.set_xticks(range(len(CATEGORIES)))
    ax.set_xticklabels(CATEGORIES, rotation=45, ha='right',
                       fontsize=PLOT_STYLE['fontsize']['tick'])
    plt.tight_layout()
    plt.savefig(output_dir / f'{sample_name}_counts.png', dpi=300)
    plt.close()

    # Size distribution plots
    for status in ['M', 'U']:
        subset = df[df['category'].str.contains(status)]
        create_distribution_plot(
            subset, 
            'Methylated' if status == 'M' else 'Unmethylated',
            sample_name,
            output_dir / f'{sample_name}_{status}_size_dist.png'
        )
        
    # Create chromosome distribution plot
    create_chromosome_distribution_plot(
        df,
        sample_name,
        output_dir / f'{sample_name}_chrom_dist.png'
    )

def process_sample(sample_name, input_dir, output_dir):
    """
    Process a single sample's methylation data
    
    Args:
        sample_name: Name of the sample to process
        input_dir: Directory containing input files
        output_dir: Directory to save outputs
    
    Returns:
        Combined DataFrame of processed data or None if error
    """
    sample_dir = Path(input_dir) / sample_name
    
    # Read data
    dfs = []
    for hap in ['hap1', 'hap2']:
        file_path = sample_dir / f"{sample_name}.{hap}.meth_regions.bed"
        df = read_methylation_file(file_path, f"Hap{hap[-1]}")
        if df is not None:
            dfs.append(df)
    
    if not dfs:
        logging.error(f"No valid data found for {sample_name}")
        return None
        
    combined_df = pd.concat(dfs)
    
    # Create plots
    create_sample_plots(combined_df, sample_name, output_dir)
    
    # Generate summary statistics
    stats = pd.DataFrame(index=CATEGORIES)
    group_stats = combined_df.groupby('category').agg({
        'size': ['count', 'mean', 'std']
    }).round(2)
    stats = group_stats.reindex(CATEGORIES, fill_value=0)
    
    # Generate chromosome-specific statistics
    chrom_stats = combined_df.groupby(['category', 'chrom']).size().unstack(fill_value=0)
    chrom_stats = chrom_stats.reindex(columns=CHROM_ORDER, fill_value=0)
    
    # Save statistics
    stats.to_csv(output_dir / f'{sample_name}_stats.txt')
    chrom_stats.to_csv(output_dir / f'{sample_name}_chrom_stats.txt')
    
    logging.info(f"Sample processing complete: {sample_name}")
    return combined_df

def main():
    import sys
    if len(sys.argv) != 4:
        print("Usage: python per_sample_analysis.py <input_dir> <output_dir> <sample_name>")
        sys.exit(1)
    
    input_dir, output_dir, sample_name = sys.argv[1:4]
    
    try:
        output_dir = Path(output_dir) 
        output_dir.mkdir(parents=True, exist_ok=True)
        process_sample(sample_name, input_dir, output_dir)
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()