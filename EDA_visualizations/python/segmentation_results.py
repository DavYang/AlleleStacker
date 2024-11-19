import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

# Check if correct number of arguments is provided
if len(sys.argv) != 4:
    print("Usage: python segmentation_results.py <input_dir> <output_dir> <sample_name>")
    sys.exit(1)

# Get input and output directories and sample name from command line arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]
sample_name = sys.argv[3]

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Define the order of labels for the x-axis
label_order = ['H1_M', 'H1_U', 'H2_M', 'H2_U', 'ASM']

# Define colors for each label
color_map = {
    'H1_M': '#800080',  # Purple
    'H1_U': '#9370DB',  # Medium purple
    'H2_M': '#006400',  # Dark green
    'H2_U': '#008000',  # Green
    'ASM': '#4B0082'    # Indigo (mix of green and blue)
}

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
    
    # Format each tick label
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

# Function to process a single file and create a plot
def process_file(file_path):
    # Read the BED file
    df = pd.read_csv(file_path, sep='\t', comment='#', header=0, names=['chrom', 'start', 'end', 'summary_label'])
    
    # Count the occurrences of each label
    counts = df['summary_label'].value_counts()
    
    # Ensure all labels are present, fill with 0 if missing
    for label in label_order:
        if label not in counts:
            counts[label] = 0
    
    # Sort the counts according to the specified order
    counts = counts.reindex(label_order)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(x=counts.index, y=counts.values, hue=counts.index, palette=color_map, legend=False, ax=ax)
    plt.title(f'Segmentation Regions for Sample: {sample_name}', fontsize=TITLE_FONT_SIZE)
    plt.xlabel('Region Label', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.ylabel('Count', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.xticks(rotation=45, fontsize=TICK_LABEL_FONT_SIZE)
    plt.yticks(fontsize=TICK_LABEL_FONT_SIZE)
    
    # Set y-axis to use scientific notation for ticks
    set_scientific_notation(ax)
    
    # Adjust y-axis label position
    ax.yaxis.set_label_coords(-0.15, 0.5)
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, f'{sample_name}_region_label_counts.png'), bbox_inches='tight')
    plt.close()

    # Generate plot for average sizes and standard deviation
    plot_region_sizes(df, sample_name, output_dir)

def plot_region_sizes(df, sample_name, output_dir):
    # Calculate region sizes
    df['size'] = df['end'] - df['start']
    
    # Calculate average size and standard deviation for each label
    size_stats = df.groupby('summary_label')['size'].agg(['mean', 'std']).reindex(label_order)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot average sizes
    sns.barplot(x=size_stats.index, y=size_stats['mean'], hue=size_stats.index, palette=color_map, legend=False, ax=ax)
    
    # Add error bars for standard deviation
    plt.errorbar(x=range(len(size_stats)), y=size_stats['mean'], yerr=size_stats['std'], fmt='none', color='black', capsize=5)
    
    plt.title(f'Average Region Sizes for Sample: {sample_name}', fontsize=TITLE_FONT_SIZE)
    plt.xlabel('Region Label', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.ylabel('Average Size (bp)', fontsize=AXIS_LABEL_FONT_SIZE)
    plt.xticks(rotation=45, fontsize=TICK_LABEL_FONT_SIZE)
    plt.yticks(fontsize=TICK_LABEL_FONT_SIZE)
    
    # Set y-axis to use scientific notation for ticks
    set_scientific_notation(ax)
    
    # Adjust y-axis label position
    ax.yaxis.set_label_coords(-0.15, 0.5)
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, f'{sample_name}_region_sizes.png'), bbox_inches='tight')
    plt.close()

# Process the file for the given sample
file_path = os.path.join(input_dir, f"{sample_name}.meth_regions.bed")
if os.path.exists(file_path):
    process_file(file_path)
    print(f"Processed {sample_name}. Plot saved in {output_dir}")
else:
    print(f"Error: File not found for sample {sample_name}")
    sys.exit(1)
