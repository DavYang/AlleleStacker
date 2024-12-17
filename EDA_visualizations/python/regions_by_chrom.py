import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os

def plot_regions_per_chromosome(h1_file_path, h2_file_path, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the BED files
    h1_df = pd.read_csv(h1_file_path, sep='\t', header=None)
    h2_df = pd.read_csv(h2_file_path, sep='\t', header=None)
    
    # Name the columns
    bed_columns = ['chrom', 'start', 'end', 'num_unmethylated', 'num_methylated', 
                  'unmethylated_samples', 'methylated_samples']
    h1_df.columns = bed_columns
    h2_df.columns = bed_columns
    
    # Count regions per chromosome for each haplotype
    h1_counts = h1_df['chrom'].value_counts().sort_index()
    h2_counts = h2_df['chrom'].value_counts().sort_index()
    
    # Get all unique chromosomes
    all_chroms = sorted(set(h1_counts.index) | set(h2_counts.index))
    
    # Create figure and axis
    plt.figure(figsize=(15, 8))
    
    # Set bar positions
    x = np.arange(len(all_chroms))
    width = 0.35
    
    # Create bars
    h1_values = [h1_counts.get(chrom, 0) for chrom in all_chroms]
    h2_values = [h2_counts.get(chrom, 0) for chrom in all_chroms]
    
    plt.bar(x - width/2, h1_values, width, label='H1', color='#3182bd', alpha=0.7)
    plt.bar(x + width/2, h2_values, width, label='H2', color='#e6550d', alpha=0.7)
    
    # Customize plot
    plt.xlabel('Chromosome', fontsize=12)
    plt.ylabel('Number of Regions', fontsize=12)
    plt.title('Distribution of Regions Across Chromosomes by Haplotype', fontsize=14, pad=20)
    plt.xticks(x, all_chroms, rotation=45)
    plt.legend()
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on top of bars
    for i, v in enumerate(h1_values):
        plt.text(i - width/2, v, str(v), ha='center', va='bottom')
    for i, v in enumerate(h2_values):
        plt.text(i + width/2, v, str(v), ha='center', va='bottom')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(output_dir, 'regions_per_chromosome.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot regions per chromosome for both haplotypes.')
    parser.add_argument('input_dir', type=str, help='Path to the input directory containing .bed files')
    parser.add_argument('output_dir', type=str, help='Directory where the plot will be saved')
    args = parser.parse_args()
    
    h1_file_path = f"{args.input_dir}/filtered_candidate_H1.bed"
    h2_file_path = f"{args.input_dir}/filtered_candidate_H2.bed"
    
    plot_regions_per_chromosome(h1_file_path, h2_file_path, args.output_dir)