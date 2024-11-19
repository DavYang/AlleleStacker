#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import random
import argparse
from typing import List, Tuple, Dict
import logging
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MethylationVisualizer:
    def __init__(self, h1_path: str, h2_path: str, output_dir: str, n_regions: int = 50):
        """
        Initialize the MethylationVisualizer with input files and parameters.
        
        Args:
            h1_path: Path to H1 haplotype BED file
            h2_path: Path to H2 haplotype BED file
            output_dir: Directory for output files
            n_regions: Number of random regions to visualize
        """
        self.h1_path = Path(h1_path)
        self.h2_path = Path(h2_path)
        self.output_dir = Path(output_dir)
        self.n_regions = n_regions
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set style parameters
        plt.style.use('seaborn')
        self.colors = {
            'unmethylated': '#3182bd',  # Blue
            'methylated': '#e6550d'      # Orange
        }
        
    def load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Load and validate the input BED files."""
        try:
            h1_df = pd.read_csv(self.h1_path, sep='\t')
            h2_df = pd.read_csv(self.h2_path, sep='\t')
            
            # Validate required columns
            required_cols = ['chrom', 'start', 'end', 'num_unmethylated', 'num_methylated']
            for df in [h1_df, h2_df]:
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    raise ValueError(f"Missing required columns: {missing_cols}")
            
            # Add chromosome order for sorting
            chrom_order = {f'chr{i}': i for i in list(range(1, 23)) + ['X', 'Y']}
            for df in [h1_df, h2_df]:
                df['chrom_order'] = df['chrom'].map(chrom_order)
                
            # Sort by chromosome and position
            h1_df = h1_df.sort_values(['chrom_order', 'start']).reset_index(drop=True)
            h2_df = h2_df.sort_values(['chrom_order', 'start']).reset_index(drop=True)
            
            return h1_df, h2_df
            
        except Exception as e:
            logger.error(f"Error loading data: {str(e)}")
            raise
    
    def select_random_regions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Select random regions from the dataset, maintaining chromosomal order."""
        if len(df) < self.n_regions:
            logger.warning(f"Requested {self.n_regions} regions but only {len(df)} available")
            return df
            
        # Stratified sampling by chromosome to ensure good genome coverage
        selected_regions = df.groupby('chrom', group_keys=False).apply(
            lambda x: x.sample(n=max(1, int(self.n_regions * len(x)/len(df))), random_state=42)
        )
        
        # If we got too many regions, subsample randomly
        if len(selected_regions) > self.n_regions:
            selected_regions = selected_regions.sample(n=self.n_regions, random_state=42)
        
        return selected_regions.sort_values(['chrom_order', 'start']).reset_index(drop=True)
    
    def create_comprehensive_plot(self, df: pd.DataFrame, haplotype: str):
        """Create a single comprehensive plot for all selected regions."""
        selected_regions = self.select_random_regions(df)
        n_regions = len(selected_regions)
        
        # Create figure
        fig_height = max(10, n_regions * 0.3)  # Adjust height based on number of regions
        fig, ax = plt.subplots(figsize=(12, fig_height))
        
        # Plot bars for each region
        y_positions = np.arange(n_regions)
        
        # Plot unmethylated bars (left side)
        ax.barh(y_positions, 
                -selected_regions['num_unmethylated'],
                height=0.6,
                color=self.colors['unmethylated'],
                alpha=0.7,
                label='Unmethylated')
        
        # Plot methylated bars (right side)
        ax.barh(y_positions,
                selected_regions['num_methylated'],
                height=0.6,
                color=self.colors['methylated'],
                alpha=0.7,
                label='Methylated')
        
        # Customize plot
        ax.set_title(f'Methylation Patterns - Haplotype {haplotype}', pad=20)
        ax.set_xlabel('Sample Count')
        
        # Create region labels
        region_labels = [f"{row['chrom']}:{row['start']:,}-{row['end']:,}" 
                        for _, row in selected_regions.iterrows()]
        ax.set_yticks(y_positions)
        ax.set_yticklabels(region_labels, fontsize=8)
        
        # Add center line
        ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5, zorder=1)
        
        # Set symmetric x-axis
        max_count = max(selected_regions['num_methylated'].max(),
                       selected_regions['num_unmethylated'].max())
        ax.set_xlim(-max_count-0.5, max_count+0.5)
        
        # Add count labels
        for idx, row in selected_regions.iterrows():
            # Unmethylated count (left)
            if row['num_unmethylated'] > 0:
                ax.text(-row['num_unmethylated']/2, idx,
                       str(row['num_unmethylated']),
                       ha='center', va='center',
                       fontsize=8)
            # Methylated count (right)
            if row['num_methylated'] > 0:
                ax.text(row['num_methylated']/2, idx,
                       str(row['num_methylated']),
                       ha='center', va='center',
                       fontsize=8)
        
        # Add legend
        ax.legend(bbox_to_anchor=(0.5, 1.05),
                 loc='center',
                 ncol=2)
        
        # Add grid
        ax.yaxis.grid(True, linestyle='--', alpha=0.3)
        
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{haplotype}_methylation_patterns.png",
                   dpi=300,
                   bbox_inches='tight')
        plt.close()
        
        # Save region data to file
        selected_regions.to_csv(self.output_dir / f"{haplotype}_selected_regions.tsv",
                              sep='\t',
                              index=False)
    
    def run(self):
        """Run the full visualization pipeline."""
        try:
            logger.info("Loading data...")
            h1_df, h2_df = self.load_data()
            
            logger.info("Creating visualization for H1 haplotype...")
            self.create_comprehensive_plot(h1_df, "H1")
            
            logger.info("Creating visualization for H2 haplotype...")
            self.create_comprehensive_plot(h2_df, "H2")
            
            logger.info("Visualization complete!")
            
        except Exception as e:
            logger.error(f"Error during visualization: {str(e)}")
            raise

def main():
    parser = argparse.ArgumentParser(description="Visualize methylation patterns across haplotypes")
    parser.add_argument("--h1", required=True, help="Path to H1 haplotype BED file")
    parser.add_argument("--h2", required=True, help="Path to H2 haplotype BED file")
    parser.add_argument("--output", required=True, help="Output directory for plots")
    parser.add_argument("--n-regions", type=int, default=50, 
                       help="Number of random regions to visualize")
    
    args = parser.parse_args()
    
    visualizer = MethylationVisualizer(
        h1_path=args.h1,
        h2_path=args.h2,
        output_dir=args.output,
        n_regions=args.n_regions
    )
    
    visualizer.run()

if __name__ == "__main__":
    main()