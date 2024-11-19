import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

def analyze_gaps(file_path, output_prefix):
    # Read the BED file
    df = pd.read_csv(file_path, sep='\t', header=0)

    # Sort the dataframe by chromosome and start position
    df = df.sort_values(['chrom', 'start'])

    # Calculate gaps between regions
    df['gap'] = df.groupby('chrom')['start'].diff()

    # Calculate mean gap and standard deviation
    mean_gap = df['gap'].mean()
    std_gap = df['gap'].std()

    # Create histogram of gaps
    plt.figure(figsize=(10, 6))
    plt.hist(df['gap'].dropna(), bins=50, edgecolor='black')
    plt.axvline(mean_gap, color='r', linestyle='dashed', linewidth=2, label=f'Mean Gap: {mean_gap:.2f}')
    plt.axvline(mean_gap - std_gap, color='g', linestyle='dashed', linewidth=2, label=f'-1 Std Dev: {mean_gap - std_gap:.2f}')
    plt.axvline(mean_gap + std_gap, color='g', linestyle='dashed', linewidth=2, label=f'+1 Std Dev: {mean_gap + std_gap:.2f}')

    plt.xlabel('Gap Size (bp)')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of Gaps Between Regions - {os.path.basename(file_path)}')
    plt.legend()
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_gap_distribution.png')
    plt.close()

    # Create box plot of gaps by chromosome
    plt.figure(figsize=(15, 8))
    df.boxplot(column='gap', by='chrom', figsize=(15, 8))
    plt.title(f'Distribution of Gaps Between Regions by Chromosome - {os.path.basename(file_path)}')
    plt.suptitle('')
    plt.ylabel('Gap Size (bp)')
    plt.yscale('log')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_gap_distribution_by_chromosome.png')
    plt.close()

    return mean_gap, std_gap

def main(hap1_dir, hap2_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    results = []

    # Process HAP1 files
    for file in os.listdir(hap1_dir):
        if file.endswith('.bed'):
            file_path = os.path.join(hap1_dir, file)
            output_prefix = os.path.join(output_dir, f"HAP1_{os.path.splitext(file)[0]}")
            mean_gap, std_gap = analyze_gaps(file_path, output_prefix)
            results.append((file, mean_gap, std_gap, "HAP1"))

    # Process HAP2 files
    for file in os.listdir(hap2_dir):
        if file.endswith('.bed'):
            file_path = os.path.join(hap2_dir, file)
            output_prefix = os.path.join(output_dir, f"HAP2_{os.path.splitext(file)[0]}")
            mean_gap, std_gap = analyze_gaps(file_path, output_prefix)
            results.append((file, mean_gap, std_gap, "HAP2"))

    # Print summary
    print("Summary of gap analysis:")
    print("HAP\tFile\tMean Gap\tStd Dev")
    for file, mean, std, hap in results:
        print(f"{hap}\t{file}\t{mean:.2f}\t{std:.2f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze gaps between regions in BED files.")
    parser.add_argument("hap1_dir", help="Directory containing HAP1 BED files")
    parser.add_argument("hap2_dir", help="Directory containing HAP2 BED files")
    parser.add_argument("output_dir", help="Directory to save output files")
    args = parser.parse_args()

    main(args.hap1_dir, args.hap2_dir, args.output_dir)