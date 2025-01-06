import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')

# Constants
PLOT_STYLE = {
    'colors': {
        'Hap1-M': '#FF0000',
        'Hap1-U': '#0000FF',
        'Hap2-M': '#FF6666',
        'Hap2-U': '#6666FF'
    },
    'fontsize': {
        'title': 20,
        'label': 15,
        'tick': 13
    }
}
CATEGORIES = ['Hap1-M', 'Hap1-U', 'Hap2-M', 'Hap2-U']

def setup_directories(input_dir, output_dir):
    """Set up input and output directories"""
    dirs = {
        'sample': Path(output_dir) / "per_sample_analysis",
        'cross': Path(output_dir) / "cross_sample_analysis"
    }
    
    if not Path(input_dir).is_dir():
        raise ValueError(f"Input directory not found: {input_dir}")
    
    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)
    
    return dirs

def read_methylation_file(file_path, haplotype):
    """Read and process methylation data file"""
    try:
        df = pd.read_csv(file_path, 
                        sep='\t',
                        comment='#',
                        names=['chrom', 'start', 'end', 'summary_label'],
                        usecols=['start', 'end', 'summary_label'])
        
        df['size'] = df['end'] - df['start']
        df['category'] = f"{haplotype}-" + df['summary_label']
        return df[['size', 'category']]  # Keep only needed columns
        
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return None

def create_distribution_plot(data, category, sample_name, output_path):
    """Create size distribution plot for a category"""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Calculate bin edges in log space
    bins = np.logspace(np.log10(data['size'].min()), np.log10(data['size'].max()), 50)
    
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

def create_sample_plots(df, sample_name, output_dir):
    """Create all sample-specific plots"""
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

def update_cross_sample_data(df, sample_name, cross_dir):
    """Update cross-sample summary data"""
    summary = {}
    
    # Calculate statistics for each category
    for hap in ['Hap1', 'Hap2']:
        for status in ['M', 'U']:
            cat = f"{hap}-{status}"
            # Create a separate DataFrame to avoid indexing issues
            subset = df[['size', 'category']].loc[df['category'] == cat].copy()
            summary.update({
                f"{hap}_{status}_count": len(subset),
                f"{hap}_{status}_mean_size": subset['size'].mean() if len(subset) > 0 else 0,
                f"{hap}_{status}_std_size": subset['size'].std() if len(subset) > 0 else 0
            })
    
    # Update or create summary file
    summary_file = cross_dir / "all_samples_summary.csv"
    summary_df = pd.DataFrame([summary], index=[sample_name])
    
    if summary_file.exists():
        existing = pd.read_csv(summary_file, index_col=0)
        # Handle duplicates by keeping the latest version
        updated = pd.concat([existing, summary_df])
        updated = updated[~updated.index.duplicated(keep='last')]
        updated.to_csv(summary_file)
    else:
        summary_df.to_csv(summary_file)
    
    return summary_file

def create_cross_sample_plots(summary_file, output_dir):
    """Create comparison plots across all samples"""
    if not summary_file.exists():
        return
        
    df = pd.read_csv(summary_file, index_col=0)
    if len(df) <= 1:
        return
    
    # Sort index to ensure consistent sample ordering
    df = df.sort_index()
    x = np.arange(len(df.index))
    width = 0.35
    
    # Create plots for different metrics
    metrics = [
        ('methylated_counts', ['Hap1_M_count', 'Hap2_M_count'], 
         'Methylated Regions Count', ['#FF0000', '#FF6666']),
        ('unmethylated_counts', ['Hap1_U_count', 'Hap2_U_count'],
         'Unmethylated Regions Count', ['#0000FF', '#6666FF'])
    ]
    
    # Plot counts
    for fname, cols, title, colors in metrics:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        bar1 = ax.bar(x - width/2, df[cols[0]], width, label='Hap1', color=colors[0])
        bar2 = ax.bar(x + width/2, df[cols[1]], width, label='Hap2', color=colors[1])
        
        ax.set_title(title, fontsize=PLOT_STYLE['fontsize']['title'])
        ax.set_xlabel('Sample', fontsize=PLOT_STYLE['fontsize']['label'])
        ax.set_ylabel('Count', fontsize=PLOT_STYLE['fontsize']['label'])
        ax.set_xticks(x)
        ax.set_xticklabels(df.index, rotation=45, ha='right', 
                          fontsize=PLOT_STYLE['fontsize']['tick'])
        ax.legend()
        plt.grid(True, alpha=0.2)
        plt.tight_layout()
        
        plt.savefig(output_dir / f'cross_sample_{fname}.png', dpi=300)
        plt.close()
    
    # Create violin plots for sizes
    size_data = []
    for hap in ['Hap1', 'Hap2']:
        mean_col = f'{hap}_M_mean_size'
        std_col = f'{hap}_M_std_size'
        count_col = f'{hap}_M_count'
        
        # Generate points based on mean and std
        for idx, row in df.iterrows():
            size_data.append({
                'Sample': idx,
                'Haplotype': hap,
                'Size': row[mean_col],
                'Std': row[std_col],
                'Count': row[count_col]
            })
    
    size_df = pd.DataFrame(size_data)
    
    # Create violin plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot violin plot with individual points
    parts = ax.violinplot([group['Size'].values for name, group in size_df.groupby('Sample')],
                         showmeans=True, showmedians=True)
    
    # Customize violin plot
    for pc in parts['bodies']:
        pc.set_facecolor('#FF9999')
        pc.set_alpha(0.6)
    parts['cmeans'].set_color('black')
    parts['cmedians'].set_color('red')
    
    ax.set_title('Methylated Region Size Distribution Across Samples', 
                fontsize=PLOT_STYLE['fontsize']['title'])
    ax.set_xlabel('Sample', fontsize=PLOT_STYLE['fontsize']['label'])
    ax.set_ylabel('Size (bp)', fontsize=PLOT_STYLE['fontsize']['label'])
    ax.set_xticks(np.arange(1, len(df.index) + 1))
    ax.set_xticklabels(df.index, rotation=45, ha='right', 
                       fontsize=PLOT_STYLE['fontsize']['tick'])
    
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_dir / 'cross_sample_size_distribution.png', dpi=300)
    plt.close()

def process_sample(sample_name, input_dir, output_dirs):
    """Process a single sample's data"""
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
    create_sample_plots(combined_df, sample_name, output_dirs['sample'])
    
    # Generate summary statistics for all categories
    stats = pd.DataFrame(index=CATEGORIES)
    group_stats = combined_df.groupby('category').agg({
        'size': ['count', 'mean', 'std']
    }).round(2)
    stats = group_stats.reindex(CATEGORIES, fill_value=0)
    
    # Save summary
    stats.to_csv(output_dirs['sample'] / f'{sample_name}_stats.txt')
    
    logging.info(f"Sample processing complete: {sample_name}")
    return combined_df

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_dir> <output_dir> <sample_name>")
        sys.exit(1)
    
    input_dir, output_dir, sample_name = sys.argv[1:4]
    
    try:
        output_dirs = setup_directories(input_dir, output_dir)
        combined_df = process_sample(sample_name, input_dir, output_dirs)
        
        if combined_df is not None:
            # Update cross-sample data and create plots
            summary_file = update_cross_sample_data(combined_df, sample_name, output_dirs['cross'])
            create_cross_sample_plots(summary_file, output_dirs['cross'])
            logging.info(f"Analysis complete for {sample_name}")
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()