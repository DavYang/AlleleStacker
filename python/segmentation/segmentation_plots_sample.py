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
from matplotlib.collections import LineCollection
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from functools import partial

# Use Agg backend for faster non-interactive plotting
plt.switch_backend('Agg')

# Set up logging with debug level
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('methylation_analysis_debug.log')
    ]
)

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
    Read and process a methylation data file with optimized settings
    
    Args:
        file_path: Path to the methylation BED file
        haplotype: Haplotype identifier (e.g., 'Hap1' or 'Hap2')
    
    Returns:
        DataFrame with processed methylation data or None if error
    """
    try:
        # Use optimized pandas read_csv settings
        df = pd.read_csv(file_path, 
                        sep='\t',
                        comment='#',
                        names=['chrom', 'start', 'end', 'summary_label'],
                        usecols=['chrom', 'start', 'end', 'summary_label'],
                        dtype={'chrom': 'category',
                               'start': np.int32,
                               'end': np.int32,
                               'summary_label': 'category'},
                        engine='c',
                        memory_map=True)
        
        # Vectorized operations
        df['size'] = df['end'] - df['start']
        df['category'] = pd.Categorical(df['summary_label'].astype(str).map(lambda x: f"{haplotype}-{x}"))
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
    Create simplified distribution plot of regions across chromosomes
    
    Args:
        df: DataFrame containing chromosome and category information
        sample_name: Name of the sample
        output_path: Path to save the plot
    """
    logging.debug(f"Creating chromosome distribution plot for {sample_name}")
    logging.debug(f"Input DataFrame shape: {df.shape}")
    logging.debug(f"Unique chromosomes: {df['chrom'].unique()}")
    logging.debug(f"Unique categories: {df['category'].unique()}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12), height_ratios=[3, 1])
    
    # Top plot: Simplified bar chart
    logging.debug("Creating top bar plot")
    bar_width = 0.2
    
    # Simple cross-tabulation of chromosome and category counts
    counts_by_chrom = pd.crosstab(df['chrom'], df['category'])
    logging.debug(f"Counts by chromosome:\n{counts_by_chrom}")
    
    # Ensure all chromosomes are included in order
    counts_by_chrom = counts_by_chrom.reindex(CHROM_ORDER, fill_value=0)
    
    # Plot bars for each category
    x = np.arange(len(CHROM_ORDER))
    for i, category in enumerate(CATEGORIES):
        if category in counts_by_chrom.columns:
            ax1.bar(x + (i * bar_width), 
                   counts_by_chrom[category], 
                   width=bar_width,
                   label=category,
                   color=PLOT_STYLE['colors'][category],
                   alpha=0.6)
    
    # Bottom plot: Simplified genome-wide distribution
    logging.debug("Creating bottom distribution plot")
    
    try:
        # Calculate chromosome positions
        chrom_positions = {}
        current_pos = 0
        
        for chrom in CHROM_ORDER:
            chrom_data = df[df['chrom'] == chrom]
            if not chrom_data.empty:
                max_pos = chrom_data['end'].max()
                chrom_positions[chrom] = (current_pos, current_pos + max_pos)
                current_pos += max_pos
        
        logging.debug(f"Chromosome positions calculated: {chrom_positions}")
        
        # Plot regions by category
        for category in CATEGORIES:
            category_data = df[df['category'] == category]
            logging.debug(f"Plotting {len(category_data)} regions for {category}")
            
            # Vectorized position calculation
            if not category_data.empty:
                # Create mapping of chromosome to start position
                chrom_starts = pd.Series({chrom: pos[0] for chrom, pos in chrom_positions.items()})
                
                # Vectorized calculation of absolute positions
                positions = (chrom_starts[category_data['chrom']].values + 
                           ((category_data['start'] + category_data['end']) / 2).values)
                
                # Create segments array directly
                segments = np.zeros((len(positions), 2, 2))
                segments[:, 0, 0] = positions  # x1
                segments[:, 1, 0] = positions  # x2
                segments[:, 1, 1] = 1          # y2 (y1 remains 0)
                
                lc = LineCollection(segments, 
                                  colors=PLOT_STYLE['colors'][category],
                                  alpha=0.1,
                                  linewidth=0.5)
                ax2.add_collection(lc)
        
        # Add chromosome boundaries
        for _, (start, end) in chrom_positions.items():
            ax2.axvline(x=end, color='black', linestyle='--', alpha=0.5)
        
        ax2.set_xlim(0, current_pos)
        ax2.set_ylim(0, 1)
        
    except Exception as e:
        logging.error(f"Error creating bottom plot: {e}")
        logging.debug("Exception details:", exc_info=True)
    
    # Customize plots
    logging.debug("Customizing plot appearance")
    
    # Top plot customization
    ax1.set_title(f'Chromosome Distribution of Methylation Regions - {sample_name}', 
                  fontsize=PLOT_STYLE['fontsize']['title'])
    ax1.set_xlabel('Chromosome', fontsize=PLOT_STYLE['fontsize']['label'])
    ax1.set_ylabel('Number of Regions', fontsize=PLOT_STYLE['fontsize']['label'])
    
    group_center = bar_width * (len(CATEGORIES) - 1) / 2
    ax1.set_xticks(x + group_center)
    ax1.set_xticklabels([chrom.replace('chr', '') for chrom in CHROM_ORDER],
                        rotation=45,
                        ha='right',
                        fontsize=PLOT_STYLE['fontsize']['tick'])
    
    ax1.grid(True, axis='y', linestyle='--', alpha=0.3)
    ax1.legend(frameon=True, 
               framealpha=1.0, 
               loc='upper right',
               fontsize=PLOT_STYLE['fontsize']['tick'])
    ax1.set_yscale('log')
    
    # Bottom plot customization
    ax2.set_title('Genome-wide Distribution of Methylation Regions', 
                  fontsize=PLOT_STYLE['fontsize']['title'])
    ax2.set_xlabel('Genomic Position', fontsize=PLOT_STYLE['fontsize']['label'])
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    plt.tight_layout()
    
    try:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logging.debug(f"Plot saved successfully to {output_path}")
    except Exception as e:
        logging.error(f"Error saving plot: {e}")
    finally:
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

def process_sample(sample_name, input_dir, base_output_dir):
    """
    Process a single sample's methylation data using parallel processing
    
    Args:
        sample_name: Name of the sample to process
        input_dir: Directory containing input files
        base_output_dir: Base directory for outputs
    """
    # Create sample-specific output directory
    output_dir = Path(base_output_dir) / sample_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Configure matplotlib for faster rendering
    plt.rcParams['path.simplify'] = True
    plt.rcParams['path.simplify_threshold'] = 1.0
    plt.rcParams['agg.path.chunksize'] = 10000
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
        
        # Use all available CPU cores except one
        num_cores = max(1, mp.cpu_count() - 1)
        
        # Process sample with parallel execution
        with ThreadPoolExecutor(max_workers=num_cores) as executor:
            process_fn = partial(process_sample, input_dir=input_dir, base_output_dir=output_dir)
            executor.submit(process_fn, sample_name)
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
