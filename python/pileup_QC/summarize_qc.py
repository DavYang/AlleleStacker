#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import sys

class CpGPlotter:
    def __init__(self, input_dir: str, output_dir: str = None):
        """Initialize plotter with input and output directories"""
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir) if output_dir else self.input_dir / 'plots'
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.output_dir / 'plotting.log'),
                logging.StreamHandler(sys.stdout)
            ]
        )

    def collect_data(self) -> pd.DataFrame:
        """Collect data from filtering summary files"""
        summary_files = list(self.input_dir.glob("*/filtering_summary.csv"))
        
        if not summary_files:
            raise ValueError(f"No filtering summary files found in {self.input_dir}")
            
        dfs = []
        for f in summary_files:
            try:
                df = pd.read_csv(f)
                # Remove duplicates if they exist
                df = df.drop_duplicates()
                dfs.append(df)
                logging.info(f"Loaded data from {f}")
            except Exception as e:
                logging.warning(f"Error reading {f}: {e}")
                
        if not dfs:
            raise ValueError("No valid data files could be loaded")
            
        combined_df = pd.concat(dfs, ignore_index=True)
        final_df = combined_df.drop_duplicates()
        
        # Log sample count information
        sample_count = len(final_df)
        logging.info(f"Analysis complete - processed {sample_count} unique samples")
        
        return final_df

    def generate_plots(self, df: pd.DataFrame) -> None:
        """Generate plots for key metrics"""
        # Common plotting parameters
        width = 12
        height = 6
        bar_width = 0.35
        colors = {'hap1': '#2ecc71', 'hap2': '#3498db'}
        
        # Custom style settings
        plt.rcParams['axes.grid'] = True
        plt.rcParams['grid.alpha'] = 0.3
        plt.rcParams['grid.linestyle'] = '--'
        
        # Define metrics to plot - matching exact column names
        metrics = {
            'Denovo': ('hap1_Denovo_CpGs', 'hap2_Denovo_CpGs', 'De Novo CpG Sites'),
            'Destroyed': ('hap1_Destroyed', 'hap2_Destroyed', 'CpG Sites Lost by Decay'),
            'Excluded': ('hap1_Excl', 'hap2_Excl', 'Total Excluded CpG Sites'),
            'Remaining': ('hap1_Kept', 'hap2_Kept', 'CpG Sites Remaining After Filter'),
            'Phantom': ('hap1_Phantom', 'hap2_Phantom', 'Phantom CpG Sites'),
            'Preserved': ('hap1_Preserved', 'hap2_Preserved', 'Preserved CpG Sites')
        }
        
        x = np.arange(len(df))
        
        # Add sample count to plot titles
        sample_count = len(df)
        
        for metric, (hap1_col, hap2_col, title) in metrics.items():
            if hap1_col not in df.columns or hap2_col not in df.columns:
                logging.warning(f"Skipping {title} plot - required columns not found")
                continue
                
            # Create figure and axis objects
            fig, ax = plt.subplots(figsize=(width, height))
            
            # Create bars
            ax.bar(x - bar_width/2, df[hap1_col], bar_width, label='Haplotype 1', 
                  color=colors['hap1'], edgecolor='black', linewidth=0.5)
            ax.bar(x + bar_width/2, df[hap2_col], bar_width, label='Haplotype 2', 
                  color=colors['hap2'], edgecolor='black', linewidth=0.5)
            
            # Calculate overall average and standard deviation
            values = np.concatenate([df[hap1_col], df[hap2_col]])
            overall_avg = np.mean(values)
            overall_std = np.std(values)
            
            # Add average line and annotation
            ax.axhline(y=overall_avg, color='red', linestyle='--', alpha=0.8, linewidth=1.5)
            ax.text(len(df)-1 + 0.5, overall_avg, 
                   f'Mean: {overall_avg:.1e} ± {overall_std:.1e}', 
                   color='red', va='center', ha='left',
                   bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
            
            # Set labels and title (now including sample count)
            ax.set_xlabel('Sample ID', fontsize=10)
            ax.set_ylabel('Count of Sites', fontsize=10)
            ax.set_title(f'{title}\n(n={sample_count} samples)', fontsize=12, pad=20)
            
            # Set x-axis ticks
            ax.set_xticks(x)
            ax.set_xticklabels(df['Sample'], rotation=45, ha='right', fontsize=8)
            
            # Format y-axis
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
            plt.yticks(fontsize=8)
            
            # Add grid
            ax.grid(True, alpha=0.3, linestyle='--')
            
            # Add legend
            ax.legend(
                frameon=True,
                fancybox=True,
                shadow=True,
                bbox_to_anchor=(0.5, -0.35),
                loc='upper center',
                borderaxespad=0.,
                ncol=2
            )
            
            # Adjust layout
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.25, right=0.85)
            
            # Save figure
            output_file = self.output_dir / f'{metric.lower()}_sites_comparison.pdf'
            plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.1)
            plt.close()
            
            # Save raw data
            stats_df = pd.DataFrame({
                'Sample': df['Sample'],
                'Haplotype 1': df[hap1_col],
                'Haplotype 2': df[hap2_col],
                'Mean': overall_avg,
                'Std Dev': overall_std
            })
            stats_file = self.output_dir / f'{metric.lower()}_stats.tsv'
            stats_df.to_csv(stats_file, sep='\t', index=False)
            logging.info(f"Generated {metric} plots and stats")

    def run(self):
        """Run complete plotting pipeline"""
        logging.info("Starting CpG plot generation")
        df = self.collect_data()
        self.generate_plots(df)
        logging.info(f"Plot generation complete. Results in {self.output_dir}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate plots from CpG filtering results')
    parser.add_argument('--input-dir', required=True,
                       help='Directory containing filtering results')
    parser.add_argument('--output-dir',
                       help='Directory for output plots (default: input_dir/plots)')
    
    args = parser.parse_args()
    
    try:
        plotter = CpGPlotter(args.input_dir, args.output_dir)
        plotter.run()
    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()