#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logging
from typing import Dict

class MultiSampleAnalyzer:
    def __init__(self, results_dir: str):
        """Initialize analyzer with directory containing sample results"""
        self.results_dir = Path(results_dir)
        self.output_dir = self.results_dir / 'multi_sample_summary'
        self.output_dir.mkdir(exist_ok=True)
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.output_dir / 'analysis.log'),
                logging.StreamHandler()
            ]
        )

    def collect_sample_data(self) -> pd.DataFrame:
        """Collect all sample results into single DataFrame"""
        summary_files = list(self.results_dir.glob("*/filtering_summary.csv"))
        if not summary_files:
            raise ValueError("No sample summary files found")
            
        dfs = []
        for f in summary_files:
            try:
                df = pd.read_csv(f)
                dfs.append(df)
            except Exception as e:
                logging.warning(f"Error reading {f}: {e}")
                
        return pd.concat(dfs, ignore_index=True)

    def generate_improved_plots(self, df: pd.DataFrame) -> None:
        """Generate separate plots for each metric"""
        # Common plotting parameters
        width = 12
        height = 6
        bar_width = 0.35
        colors = {'hap1': '#2ecc71', 'hap2': '#3498db'}
        
        # Custom style settings
        plt.rcParams['axes.grid'] = True
        plt.rcParams['grid.alpha'] = 0.3
        plt.rcParams['grid.linestyle'] = '--'
        
        # Calculate adjusted destroyed sites
        df['hap1_AdjDestroyed'] = df['hap1_Destroyed'] - df['hap1_Preserved']
        df['hap2_AdjDestroyed'] = df['hap2_Destroyed'] - df['hap2_Preserved']
        
        metrics = {
            'Kept': ('hap1_Kept', 'hap2_Kept', 'CpG Sites Retained'),
            'Excl': ('hap1_Excl', 'hap2_Excl', 'CpG Sites Removed'),
            'Destroyed': ('hap1_Destroyed', 'hap2_Destroyed', 'CpG Sites Destroyed'),
            'Phantom': ('hap1_Phantom', 'hap2_Phantom', 'Phantom CpG Sites'),
            'Preserved': ('hap1_Preserved', 'hap2_Preserved', 'Preserved CpG Sites'),
            'AdjDestroyed': ('hap1_AdjDestroyed', 'hap2_AdjDestroyed', 'Adjusted CpG Sites Destroyed (Destroyed - Preserved)')
        }
        
        x = np.arange(len(df))
        
        for metric, (hap1_col, hap2_col, title) in metrics.items():
            if hap1_col not in df.columns or hap2_col not in df.columns:
                logging.warning(f"Skipping {title} plot - columns not found in data")
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
            
            # Add average line and annotation with standard deviation
            ax.axhline(y=overall_avg, color='red', linestyle='--', alpha=0.8, linewidth=1.5)
            ax.text(len(df)-1 + 0.5, overall_avg, 
                   f'Group Mean: {overall_avg:.1e} ± {overall_std:.1e}', 
                   color='red', va='center', ha='left',
                   bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
            
            # Set labels and title
            ax.set_xlabel('Sample ID', fontsize=10)
            ax.set_ylabel('Count of Sites', fontsize=10)
            ax.set_title(title, fontsize=12, pad=20)
            
            # Set x-axis ticks
            ax.set_xticks(x)
            ax.set_xticklabels(df['Sample'], rotation=45, ha='right', fontsize=8)
            
            # Format y-axis with scientific notation
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1e}'))
            plt.yticks(fontsize=8)
            
            # Add grid
            ax.grid(True, alpha=0.3, linestyle='--')
            
            # Add legend below x-axis label
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
            plt.savefig(
                self.output_dir / f'{metric.lower()}_sites_comparison.pdf',
                dpi=300,
                bbox_inches='tight',
                pad_inches=0.1
            )
            plt.close()

    def calculate_statistics(self, df: pd.DataFrame) -> pd.DataFrame:
        """Calculate summary statistics for each metric"""
        # Calculate adjusted destroyed sites
        df['hap1_AdjDestroyed'] = df['hap1_Destroyed'] - df['hap1_Preserved']
        df['hap2_AdjDestroyed'] = df['hap2_Destroyed'] - df['hap2_Preserved']
        
        stats_data = []
        metrics = ['Total', 'Kept', 'Excl', 'Destroyed', 'Phantom', 'Preserved', 'AdjDestroyed']
        
        for hap in ['hap1', 'hap2']:
            for metric in metrics:
                col = f'{hap}_{metric}'
                if col in df.columns:
                    stats = df[col].describe()
                    stats_data.append({
                        'Metric': metric,
                        'Haplotype': hap,
                        'Mean': stats['mean'],
                        'Std': stats['std'],
                        'Min': stats['min'],
                        'Max': stats['max'],
                        'Median': stats['50%'],
                        'Count': stats['count']
                    })
        
        return pd.DataFrame(stats_data)

    def write_summary_report(self, df: pd.DataFrame, stats_df: pd.DataFrame) -> None:
        """Generate detailed summary report"""
        report_path = self.output_dir / 'summary_report.txt'
        
        with open(report_path, 'w') as f:
            f.write("CpG Analysis Summary Report\n")
            f.write("==========================\n\n")
            
            # Write statistics for each haplotype
            for hap in ['hap1', 'hap2']:
                f.write(f"\n{hap.upper()} Summary Statistics:\n")
                f.write("-" * 30 + "\n")
                
                hap_stats = stats_df[stats_df['Haplotype'] == hap]
                for _, row in hap_stats.iterrows():
                    f.write(f"\n{row['Metric']}:\n")
                    f.write(f"  Mean ± Std: {row['Mean']:,.2f} ± {row['Std']:,.2f}\n")
                    f.write(f"  Range: {row['Min']:,.2f} - {row['Max']:,.2f}\n")
                    f.write(f"  Median: {row['Median']:,.2f}\n")
            
            # Calculate and write ratios
            f.write("\nKey Ratios:\n")
            f.write("-" * 20 + "\n")
            for hap in ['hap1', 'hap2']:
                kept = df[f'{hap}_Kept'].mean()
                total = df[f'{hap}_Total'].mean()
                excl = df[f'{hap}_Excl'].mean()
                destroyed = df[f'{hap}_Destroyed'].mean()
                preserved = df[f'{hap}_Preserved'].mean()
                adj_destroyed = df[f'{hap}_AdjDestroyed'].mean()
                
                f.write(f"\n{hap.upper()}:\n")
                f.write(f"  Kept/Total: {(kept/total)*100:.2f}%\n")
                f.write(f"  Excluded/Total: {(excl/total)*100:.2f}%\n")
                f.write(f"  Adjusted Destroyed/Total: {(adj_destroyed/total)*100:.2f}%\n")

    def run_analysis(self) -> None:
        """Run complete analysis pipeline"""
        logging.info("Starting multi-sample analysis")
        
        # Collect and process data
        df = self.collect_sample_data()
        stats_df = self.calculate_statistics(df)
        
        # Save processed data
        df.to_csv(self.output_dir / 'all_samples_summary.tsv', sep='\t', index=False)
        stats_df.to_csv(self.output_dir / 'aggregate_statistics.tsv', sep='\t', index=False)
        
        # Generate visualizations
        self.generate_improved_plots(df)
        
        # Write summary report
        self.write_summary_report(df, stats_df)
        
        logging.info(f"Analysis complete. Results in {self.output_dir}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate multi-sample summary for CpG filtering results')
    parser.add_argument('--results-dir', required=True, 
                       help='Directory containing individual sample results')
    
    args = parser.parse_args()
    
    try:
        analyzer = MultiSampleAnalyzer(args.results_dir)
        analyzer.run_analysis()
    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")
        raise

if __name__ == '__main__':
    main()