import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse

def plot_distribution(input_dir, output_dir, sample_name):
    h1_file = os.path.join(input_dir, "regions", "H1_U", f"{sample_name}_H1_U.bed")
    h2_file = os.path.join(input_dir, "regions", "H2_U", f"{sample_name}_H2_U.bed")
    dfs = []

    for file, hap in [(h1_file, "H1"), (h2_file, "H2")]:
        if os.path.exists(file):
            try:
                df = pd.read_csv(file, sep='\t')
                df['haplotype'] = hap
                dfs.append(df)
            except pd.errors.EmptyDataError:
                print(f"Warning: {file} is empty. Skipping.")
            except pd.errors.ParserError as e:
                print(f"Error reading {file}: {e}")

    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        plt.figure(figsize=(12, 6))
        sns.histplot(data=combined_df, x='size', hue='haplotype', kde=True, bins=50, log_scale=True)
        plt.title(f"Distribution of Unmethylated Region Sizes for {sample_name}")
        plt.xlabel("Region Size (bp)")
        plt.ylabel("Frequency")
        output_path = os.path.join(output_dir, f"{sample_name}_unmethylated_region_distribution.png")
        plt.savefig(output_path)
        plt.close()
        print(f"Plot saved to: {output_path}")
    else:
        print(f"No unmethylated regions found for sample {sample_name}")


def main():
    parser = argparse.ArgumentParser(description="Plot the distribution of unmethylated regions.")
    parser.add_argument("input_dir", help="Path to the input directory containing the extracted regions.")
    parser.add_argument("output_dir", help="Path to the output directory for the plots.")
    parser.add_argument("sample_name", help="Name of the sample to process.")
    args = parser.parse_args()

    plot_distribution(args.input_dir, args.output_dir, args.sample_name)

if __name__ == "__main__":
    main()
