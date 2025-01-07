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

def process_sample(sample_name, input_dir, output_dir):
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
    create_sample_plots(combined_df, sample_name, output_dir)
    
    # Generate summary statistics for all categories
    stats = pd.DataFrame(index=CATEGORIES)
    group_stats = combined_df.groupby('category').agg({
        'size': ['count', 'mean', 'std']
    }).round(2)
    stats = group_stats.reindex(CATEGORIES, fill_value=0)
    
    # Save summary
    stats.to_csv(output_dir / f'{sample_name}_stats.txt')
    
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