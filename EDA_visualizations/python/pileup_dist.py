import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker

def process_sample(sample_dir, output_dir):
    sample_name = os.path.basename(os.path.normpath(sample_dir))
    
    def read_file_in_chunks(file_path, chunksize=1000000):
        for chunk in pd.read_csv(file_path, sep='\t', header=None, 
                                 names=['chr', 'start', 'end', 'meth_score', 'type', 'total', 'meth', 'unmeth', 'meth_percent'],
                                 usecols=['meth_score', 'meth_percent'], chunksize=chunksize):
            yield chunk

    summary_stats = {}

    for file_suffix in ['combined', 'hap1', 'hap2']:
        file_path = os.path.join(sample_dir, f"{sample_name}.GRCh38.{file_suffix}.bed")
        if os.path.exists(file_path):
            print(f"Processing file: {file_path}")
            meth_scores = []
            meth_percentages = []
            for chunk in read_file_in_chunks(file_path):
                meth_scores.extend(chunk['meth_score'])
                meth_percentages.extend(chunk['meth_percent'])

            meth_scores = np.array(meth_scores)
            meth_percentages = np.array(meth_percentages)

            # Create plot with reduced size
            plt.figure(figsize=(10, 6))
            
            # Increase font size
            plt.rcParams.update({'font.size': 12})
            
            # Plot distribution of methylation scores
            plt.hist(meth_scores, bins=range(0, 105, 5), edgecolor='black')
            plt.title(f'Methylation Probability Pileup Distribution: {sample_name} ({file_suffix})', fontsize=20)
            plt.xlabel('Methylation Score', fontsize=15)
            plt.ylabel('# of CpGs', fontsize=15)
            plt.xlim(0, 100)

            # Set y-axis to scientific notation on tick labels
            ax = plt.gca()
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{x:.1e}'))

            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{sample_name}_{file_suffix}_methylation_scores_distribution.png'), dpi=300)
            plt.close()

            
            # Generate summary statistics
            summary_stats[file_suffix] = {
                "Methylation Scores": pd.Series(meth_scores).describe(),
                "Methylation Percentages": pd.Series(meth_percentages).describe()
            }
        else:
            print(f"Warning: File not found: {file_path}")

    return summary_stats

def main(sample_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Sample directory: {sample_dir}")
    print(f"Output directory: {output_dir}")
    print("Contents of sample directory:")
    for item in os.listdir(sample_dir):
        print(f"- {item}")

    print(f"Processing {os.path.basename(sample_dir)}...")
    summary_stats = process_sample(sample_dir, output_dir)

    # Save summary statistics to a file
    with open(os.path.join(output_dir, f'{os.path.basename(sample_dir)}_summary_statistics.txt'), 'w') as f:
        f.write(f"Summary Statistics for {os.path.basename(sample_dir)}:\n")
        for file_suffix, stats in summary_stats.items():
            f.write(f"\n{file_suffix.upper()}:\n")
            for stat_type, stat_values in stats.items():
                f.write(f"\n{stat_type}:\n")
                f.write(str(stat_values))
            f.write("\n" + "="*50 + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pileup_dist.py <sample_directory> <output_directory>")
        sys.exit(1)
    
    sample_dir = sys.argv[1]
    output_dir = sys.argv[2]
    main(sample_dir, output_dir)
