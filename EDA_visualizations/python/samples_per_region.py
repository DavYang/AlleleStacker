import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import logging

def setup_logger():
    """Configure logging"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def add_value_labels(ax, rects, y_offset=0.1):
    """Add value labels with vertical orientation"""
    for rect in rects:
        height = rect.get_height()
        if height > 0:  # Only add label if bar height > 0
            ax.text(rect.get_x() + rect.get_width()/2.,
                   height + (height * y_offset),
                   f'{int(height):,}',
                   ha='center',
                   va='bottom',
                   fontsize=8,
                   rotation=90)

def create_methylation_plots(df_h1, df_h2, output_prefix, logger):
    """Create and save plots showing distribution of region counts"""
    logger.info("Creating methylation distribution plots...")
    
    # Set up the figure
    plt.rcParams['figure.figsize'] = [15, 6]
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 11
    
    # Create figure and axes
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.patch.set_facecolor('white')

    # Get value counts for unmethylated and methylated samples
    h1_unmeth_counts = df_h1['num_unmethylated'].value_counts().sort_index()
    h2_unmeth_counts = df_h2['num_unmethylated'].value_counts().sort_index()
    h1_meth_counts = df_h1['num_methylated'].value_counts().sort_index()
    h2_meth_counts = df_h2['num_methylated'].value_counts().sort_index()
    
    # Set bar properties
    bar_width = 0.35
    
    # Plot 1: Unmethylated samples distribution
    max_unmeth = max(max(h1_unmeth_counts.index), max(h2_unmeth_counts.index))
    x = np.arange(max_unmeth + 1)
    
    rects1 = ax1.bar(x - bar_width/2, [h1_unmeth_counts.get(i, 0) for i in x], 
                     bar_width, label='H1', color='#3498db', alpha=0.7)
    rects2 = ax1.bar(x + bar_width/2, [h2_unmeth_counts.get(i, 0) for i in x],
                     bar_width, label='H2', color='#2ecc71', alpha=0.7)
    
    # Add labels to unmethylated plot
    for rect in rects1:
        add_value_labels(ax1, [rect], y_offset=0.05)
    for rect in rects2:
        add_value_labels(ax1, [rect], y_offset=0.05)
    
    ax1.set_xlabel('Number of Unmethylated Samples per Region')
    ax1.set_ylabel('Count of Regions')
    ax1.set_title('Distribution of Regions by\nNumber of Unmethylated Samples')
    ax1.legend(frameon=False)
    
    # Set axis properties for unmethylated plot
    ax1.set_xticks(x)
    ax1.set_xlim(-0.5, max_unmeth + 0.5)
    ax1.set_ylim(0, 5000)
    ax1.set_yticks(np.arange(0, 5001, 500))
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.tick_params(direction='out')
    
    # Plot 2: Methylated samples distribution
    max_meth = max(max(h1_meth_counts.index), max(h2_meth_counts.index))
    x = np.arange(max_meth + 1)
    
    rects3 = ax2.bar(x - bar_width/2, [h1_meth_counts.get(i, 0) for i in x],
                     bar_width, label='H1', color='#3498db', alpha=0.7)
    rects4 = ax2.bar(x + bar_width/2, [h2_meth_counts.get(i, 0) for i in x],
                     bar_width, label='H2', color='#2ecc71', alpha=0.7)
    
    # Add labels to methylated plot
    for rect in rects3:
        add_value_labels(ax2, [rect], y_offset=0.05)
    for rect in rects4:
        add_value_labels(ax2, [rect], y_offset=0.05)
    
    ax2.set_xlabel('Number of Methylated Samples per Region')
    ax2.set_ylabel('Count of Regions')
    ax2.set_title('Distribution of Regions by\nNumber of Methylated Samples')
    ax2.legend(frameon=False)
    
    # Set axis properties for methylated plot
    ax2.set_xticks(x)
    ax2.set_xlim(-0.5, max_meth + 0.5)
    ax2.set_ylim(0, 3500)
    ax2.set_yticks(np.arange(0, 3501, 500))
    
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(direction='out')

    plt.tight_layout()
    
    # Save plot
    os.makedirs(os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else '.', exist_ok=True)
    plot_path = f'{output_prefix}_histograms.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    logger.info(f"Saved plot to {plot_path}")
    return plot_path

def save_statistics(df_h1, df_h2, output_prefix, logger):
    """Save summary statistics"""
    logger.info("Generating statistics summary...")
    
    stats_path = f'{output_prefix}_statistics.txt'
    with open(stats_path, 'w') as f:
        f.write("=== Methylation Distribution Statistics ===\n\n")
        
        # Write counts for each haplotype
        for name, df in [('H1', df_h1), ('H2', df_h2)]:
            f.write(f"\n=== {name} Statistics ===\n")
            f.write(f"Total regions: {len(df):,}\n")
            f.write(f"Mean unmethylated samples per region: {df['num_unmethylated'].mean():.2f}\n")
            f.write(f"Mean methylated samples per region: {df['num_methylated'].mean():.2f}\n")
            f.write(f"Median unmethylated samples per region: {df['num_unmethylated'].median():.2f}\n")
            f.write(f"Median methylated samples per region: {df['num_methylated'].median():.2f}\n")
    
    logger.info(f"Saved statistics to {stats_path}")
    return stats_path

def main():
    parser = argparse.ArgumentParser(description='Analyze methylation distribution patterns')
    parser.add_argument('--h1-input', required=True, help='Path to H1 input file')
    parser.add_argument('--h2-input', required=True, help='Path to H2 input file')
    parser.add_argument('--output-prefix', required=True, help='Prefix for output files')
    
    args = parser.parse_args()
    logger = setup_logger()
    
    try:
        logger.info("Reading input files...")
        df_h1 = pd.read_csv(args.h1_input, sep='\t')
        df_h2 = pd.read_csv(args.h2_input, sep='\t')
        
        plot_path = create_methylation_plots(df_h1, df_h2, args.output_prefix, logger)
        stats_path = save_statistics(df_h1, df_h2, args.output_prefix, logger)
        
        logger.info("Analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Error occurred: {str(e)}")
        raise

if __name__ == "__main__":
    main()