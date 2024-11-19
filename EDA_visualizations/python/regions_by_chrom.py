import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np

def plot_regions_per_chromosome(h1_file_path, h2_file_path, output_prefix):
    # Read the data for both haplotypes
    df_h1 = pd.read_csv(h1_file_path, sep='\t')
    df_h2 = pd.read_csv(h2_file_path, sep='\t')
    
    # Count the number of regions per chromosome for both haplotypes
    chromosome_counts_h1 = df_h1['chrom'].value_counts().sort_index()
    chromosome_counts_h2 = df_h2['chrom'].value_counts().sort_index()
    
    # Calculate overall averages
    h1_avg = chromosome_counts_h1.mean()
    h2_avg = chromosome_counts_h2.mean()
    
    # Define the order of chromosomes
    chromosome_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    
    # Reindex the counts to ensure all chromosomes are included and in the correct order
    chromosome_counts_h1 = chromosome_counts_h1.reindex(chromosome_order, fill_value=0)
    chromosome_counts_h2 = chromosome_counts_h2.reindex(chromosome_order, fill_value=0)
    
    # Calculate per-chromosome averages
    chr_averages = (chromosome_counts_h1 + chromosome_counts_h2) / 2
    
    # Set up the plot with larger figure size
    plt.figure(figsize=(15, 8))
    
    # Set font sizes
    SMALL_SIZE = 15
    MEDIUM_SIZE = 18
    LARGE_SIZE = 26

    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=LARGE_SIZE)
    plt.rc('axes', labelsize=MEDIUM_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=MEDIUM_SIZE)
    
    # Set the positions of the bars
    x = np.arange(len(chromosome_order))
    width = 0.35
    
    # Create the bars
    plt.bar(x - width/2, chromosome_counts_h1, width, 
            label=f'Haplotype 1', color='purple', alpha=0.7)
    plt.bar(x + width/2, chromosome_counts_h2, width, 
            label=f'Haplotype 2', color='green', alpha=0.7)
    
    # Add average values as diagonal text above the bars
    for i, avg in enumerate(chr_averages):
        plt.text(x[i], max(chromosome_counts_h1[i], chromosome_counts_h2[i]) + 50,
                f'{avg:.1f}', ha='center', va='bottom', fontsize=12,
                rotation=30)  # Added rotation for diagonal text
    
    # Customize the plot
    # plt.xlabel('Chromosome', fontsize=MEDIUM_SIZE, labelpad=10)
    plt.ylabel('# of Regions', fontsize=MEDIUM_SIZE, labelpad=10)
    plt.title('Genome Wide Distribution of Allele Stacker Consensus Regions', fontsize=LARGE_SIZE, pad=20)
    plt.xticks(x, chromosome_order, rotation=45)
    plt.legend()
    
    # Remove top and right spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    # Add grid lines for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Set y-axis limit without extra top tick
    max_count = max(max(chromosome_counts_h1), max(chromosome_counts_h2))
    y_max = max_count * 1.15  # Increased padding to accommodate average labels
    plt.ylim(0, y_max)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plot_path = f'{output_prefix}_regions_per_chromosome.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Haplotype 1 average regions per chromosome: {h1_avg:.1f}")
    print(f"Haplotype 2 average regions per chromosome: {h2_avg:.1f}")
    print("\nPer-chromosome averages:")
    for chr_name, avg in zip(chromosome_order, chr_averages):
        print(f"{chr_name}: {avg:.1f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot regions per chromosome for both haplotypes.')
    parser.add_argument('input_dir', type=str, help='Path to the input directory containing .bed files')
    parser.add_argument('output_prefix', type=str, help='Prefix for the output plot file')
    args = parser.parse_args()
    
    h1_file_path = f"{args.input_dir}/filtered_consensus_H1.bed"
    h2_file_path = f"{args.input_dir}/filtered_consensus_H2.bed"
    
    plot_regions_per_chromosome(h1_file_path, h2_file_path, args.output_prefix)