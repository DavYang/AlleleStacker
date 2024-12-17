import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

# Check if correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: python methylation_analysis.py <sample_name>")
    print("Example: python methylation_analysis.py SPM180")
    sys.exit(1)

# Get sample name from command line arguments
sample_name = sys.argv[1]

# Files are in current directory
current_dir = os.getcwd()

# Define the names of input files based on your structure
hap1_file = f"{sample_name}.hap1.meth_regions.bed"
hap2_file = f"{sample_name}.hap2.meth_regions.bed"

# Define colors for visualization - red for methylated, blue for unmethylated
color_map = {
    'Hap1-M': '#FF0000',  # Red
    'Hap1-U': '#0000FF',  # Blue
    'Hap2-M': '#FF6666',  # Light red
    'Hap2-U': '#6666FF'   # Light blue
}

# Define the order of categories for consistent plotting
category_order = ['Hap1-M', 'Hap1-U', 'Hap2-M', 'Hap2-U']

# Define font sizes
TITLE_FONT_SIZE = 20
AXIS_LABEL_FONT_SIZE = 15
TICK_LABEL_FONT_SIZE = 13

def set_scientific_notation(ax):
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

def read_and_process_file(file_path, haplotype):
    """Read and process a single methylation regions file"""
    if not os.path.exists(file_path):
        print(f"Error: File not found: {file_path}")
        return None
    
    df = pd.read_csv(file_path, sep='\t', comment='#')
    df['size'] = df['end'] - df['start']
    # Add haplotype-specific labels
    df['category'] = df['summary_label'].apply(lambda x: f"{haplotype}-{x}")
    return df

def create_plots(combined_df, sample_name):
    """Create plots for methylation status across haplotypes"""
    
    # Plot 1: Region Counts by Methylation Status and Haplotype
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Count regions for each category
    counts = combined_df['category'].value_counts().reindex(category_order).fillna(0)
    
    sns.barplot(x=counts.index, y=counts.values, palette=[color_map[cat] for cat in category_order], ax=ax)
    plt.title(f'Methylation Region Counts for {sample_name}', fontsize=TITLE_FONT_SIZE)
    plt.xlabel('Category', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.ylabel('Count', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.xticks(rotation=45, fontsize=TICK_LABEL_FONT_SIZE)
    plt.yticks(fontsize=TICK_LABEL_FONT_SIZE)
    
    set_scientific_notation(ax)
    ax.yaxis.set_label_coords(-0.1, 0.5)
    
    plt.savefig(f'{sample_name}_methylation_counts.png', bbox_inches='tight')
    plt.close()
    
    # Plot 2: Region Sizes by Methylation Status and Haplotype
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Calculate statistics for each category
    size_stats = combined_df.groupby('category')['size'].agg(['mean', 'std']).reindex(category_order)
    
    # Create bar plot with error bars
    sns.barplot(x=size_stats.index, y=size_stats['mean'], 
               palette=[color_map[cat] for cat in category_order], ax=ax)
    plt.errorbar(x=range(len(size_stats)), y=size_stats['mean'], 
                yerr=size_stats['std'], fmt='none', color='black', capsize=5)
    
    plt.title(f'Average Region Sizes for {sample_name}', fontsize=TITLE_FONT_SIZE)
    plt.xlabel('Category', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.ylabel('Average Size (bp)', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.xticks(rotation=45, fontsize=TICK_LABEL_FONT_SIZE)
    plt.yticks(fontsize=TICK_LABEL_FONT_SIZE)
    
    set_scientific_notation(ax)
    ax.yaxis.set_label_coords(-0.1, 0.5)
    
    plt.savefig(f'{sample_name}_region_sizes.png', bbox_inches='tight')
    plt.close()

def print_summary_statistics(combined_df):
    """Print summary statistics for all categories"""
    print("\nSummary Statistics:")
    
    # Count statistics
    counts = combined_df['category'].value_counts().reindex(category_order).fillna(0)
    print("\nRegion Counts:")
    for category in category_order:
        print(f"{category}: {counts[category]:.0f} regions")
    
    # Size statistics
    size_stats = combined_df.groupby('category')['size'].agg(['mean', 'std']).reindex(category_order)
    print("\nRegion Sizes:")
    for category in category_order:
        if category in size_stats.index:
            mean = size_stats.loc[category, 'mean']
            std = size_stats.loc[category, 'std']
            print(f"{category}: {mean:.2f} ± {std:.2f} bp")

# Main execution
def main():
    # Read and process both files
    hap1_df = read_and_process_file(hap1_file, "Hap1")
    hap2_df = read_and_process_file(hap2_file, "Hap2")
    
    if hap1_df is None or hap2_df is None:
        sys.exit(1)
    
    # Combine the dataframes
    combined_df = pd.concat([hap1_df, hap2_df])
    
    # Create the plots
    create_plots(combined_df, sample_name)
    print(f"Analysis complete for {sample_name}. Plots saved in current directory.")
    
    # Print summary statistics
    print_summary_statistics(combined_df)

if __name__ == "__main__":
    main()