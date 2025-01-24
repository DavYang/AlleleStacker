import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path

def plot_distribution(input_dir, output_dir, sample_name):
    """
    Create distribution plots for methylated and unmethylated regions.
    
    Args:
        input_dir: Directory containing the H1_M, H1_U, H2_M, H2_U folders
        output_dir: Directory to save the plots
        sample_name: Name of the sample to process
    """
    # Initialize data collection
    dfs = []
    
    # Process each methylation state
    for hap in ['H1', 'H2']:
        for state in ['M', 'U']:
            file_path = Path(input_dir) / f"{hap}_{state}" / f"{sample_name}_{hap}_{state}.bed"
            if file_path.exists():
                try:
                    df = pd.read_csv(file_path, sep='\t')
                    df['size'] = df['end'] - df['start']
                    df['haplotype'] = hap
                    df['state'] = state
                    dfs.append(df)
                except pd.errors.EmptyDataError:
                    print(f"Warning: {file_path} is empty")
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
    
    if not dfs:
        print(f"No data found for sample {sample_name}")
        return
    
    # Combine all data
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Create separate plots for methylated and unmethylated regions
    for state in ['M', 'U']:
        state_df = combined_df[combined_df['state'] == state]
        if not state_df.empty:
            plt.figure(figsize=(12, 6))
            sns.histplot(data=state_df, 
                        x='size', 
                        hue='haplotype',
                        multiple="layer",
                        palette={'H1': '#FF0000' if state == 'M' else '#0000FF',
                                'H2': '#FF6666' if state == 'M' else '#6666FF'},
                        alpha=0.5,
                        log_scale=True,
                        bins=50)
            
            plt.title(f"Distribution of {'Methylated' if state == 'M' else 'Unmethylated'} "
                     f"Region Sizes - {sample_name}")
            plt.xlabel("Region Size (bp)")
            plt.ylabel("Count")
            
            # Add grid
            plt.grid(True, which="major", ls="-", alpha=0.2)
            plt.grid(True, which="minor", ls=":", alpha=0.1)
            
            # Improve legend
            plt.legend(title="Haplotype", frameon=True, framealpha=1.0)
            
            # Save plot
            output_path = Path(output_dir) / f"{sample_name}_{state}_region_distribution.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot distribution of methylated/unmethylated regions.")
    parser.add_argument("--input-dir", required=True, help="Path to the input directory containing region folders")
    parser.add_argument("--output-dir", required=True, help="Path to the output directory for plots")
    parser.add_argument("--sample-name", required=True, help="Name of the sample to process")
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create plots
    plot_distribution(args.input_dir, args.output_dir, args.sample_name)

if __name__ == "__main__":
    main()
