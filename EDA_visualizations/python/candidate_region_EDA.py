import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Analyze BED files for region distributions')
    parser.add_argument('--input-dir', type=str, required=True, help='Directory containing input BED files')
    parser.add_argument('--output-dir', type=str, required=True, help='Directory for output plots')
    return parser.parse_args()

def load_bed_file(file_path):
    """Load BED file and return DataFrame."""
    return pd.read_csv(file_path, sep='\t')

def plot_sample_count_distribution(df, output_path, haplotype):
    """Plot how many regions have specific numbers of samples."""
    # Calculate total samples per region (methylated + unmethylated)
    df['total_samples'] = df['num_methylated'] + df['num_unmethylated']
    
    # Count frequency of each total
    sample_counts = df['total_samples'].value_counts().sort_index()
    
    plt.figure(figsize=(12, 6))
    plt.bar(sample_counts.index, sample_counts.values, color='skyblue')
    plt.xlabel('Number of Samples per Region')
    plt.ylabel('Number of Regions')
    plt.title(f'Distribution of Samples per Region (H{haplotype})')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add count labels on top of bars
    for i, v in zip(sample_counts.index, sample_counts.values):
        plt.text(i, v, str(v), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(output_path / f'sample_distribution_H{haplotype}.png')
    plt.close()

def plot_region_size_distribution(df, output_path, haplotype):
    """Plot distribution of region sizes."""
    plt.figure(figsize=(12, 6))
    
    # Plot histogram with log scale for better visualization
    plt.hist(df['size'], bins=50, color='lightgreen', edgecolor='black')
    plt.xlabel('Region Size (bp)')
    plt.ylabel('Count')
    plt.title(f'Distribution of Region Sizes (H{haplotype})')
    plt.yscale('log')  # Log scale for y-axis due to likely skewed distribution
    
    # Add summary statistics as text
    stats = df['size'].describe()
    stats_text = f'Mean: {stats["mean"]:.0f}\nMedian: {stats["50%"]:.0f}\nMax: {stats["max"]:.0f}'
    plt.text(0.95, 0.95, stats_text,
             transform=plt.gca().transAxes,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path / f'region_size_distribution_H{haplotype}.png')
    plt.close()

def plot_chromosome_distribution(df, output_path, haplotype):
    """Plot distribution of regions across chromosomes (autosomes + sex chromosomes)."""
    # Create ordered list of chromosomes
    autosomes = [f'chr{i}' for i in range(1, 23)]
    sex_chromosomes = ['chrX', 'chrY']
    ordered_chromosomes = autosomes + sex_chromosomes
    
    # Count regions per chromosome
    chrom_counts = df['chrom'].value_counts()
    
    # Create ordered data for plotting
    ordered_counts = []
    ordered_labels = []
    for chrom in ordered_chromosomes:
        if chrom in chrom_counts.index:
            ordered_counts.append(chrom_counts[chrom])
            # Remove 'chr' prefix for x-axis labels
            ordered_labels.append(chrom.replace('chr', ''))
    
    plt.figure(figsize=(15, 6))
    bars = plt.bar(range(len(ordered_counts)), ordered_counts, color='salmon')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of Regions')
    plt.title(f'Distribution of Regions Across Chromosomes (H{haplotype})')
    
    # Set x-axis labels
    plt.xticks(range(len(ordered_counts)), ordered_labels)
    
    # Add count labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height):,}',
                ha='center', va='bottom')
    
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path / f'chromosome_distribution_H{haplotype}.png')
    plt.close()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Create Path objects
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each haplotype
    for haplotype in ['1', '2']:
        print(f"\nProcessing Haplotype {haplotype}...")
        
        # Load data
        file_path = input_dir / f'filtered_candidate_H{haplotype}.bed'
        df = load_bed_file(file_path)
        
        # Generate plots
        plot_sample_count_distribution(df, output_dir, haplotype)
        plot_region_size_distribution(df, output_dir, haplotype)
        plot_chromosome_distribution(df, output_dir, haplotype)
        
        # Print summary statistics
        print(f"Total regions: {len(df):,}")
        print(f"Average region size: {df['size'].mean():.1f} bp")
        print(f"Median region size: {df['size'].median():.1f} bp")
        print(f"Total bases covered: {df['size'].sum():,} bp")

if __name__ == "__main__":
    main()