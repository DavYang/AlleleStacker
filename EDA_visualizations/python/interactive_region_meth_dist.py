#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import plotly.graph_objects as go
import plotly.subplots as sp
import argparse
from typing import Tuple, Dict
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
        """Initialize the MethylationVisualizer with input files and parameters."""
        self.h1_path = Path(h1_path)
        self.h2_path = Path(h2_path)
        self.output_dir = Path(output_dir)
        self.n_regions = n_regions
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.colors = {
            'unmethylated': '#3182bd',  # Blue
            'methylated': '#e6550d'      # Orange
        }
        
        # Figure dimensions
        self.fig_width = 1200
        self.fig_height_per_region = 25  # Height per region
        self.min_fig_height = 800        # Minimum figure height
    
    def calculate_methylation_metrics(self, row: pd.Series) -> Dict[str, float]:
        """Calculate methylation metrics for a single region."""
        total_samples = row['num_methylated'] + row['num_unmethylated']
        
        if total_samples == 0:
            return {
                'total_samples': 0,
                'meth_proportion': 0.0,
                'unmeth_proportion': 0.0,
                'proportion_diff': 0.0,
                'sample_weight': 0.0,
                'final_score': 0.0
            }
            
        meth_proportion = row['num_methylated'] / total_samples
        unmeth_proportion = row['num_unmethylated'] / total_samples
        proportion_diff = meth_proportion - unmeth_proportion
        sample_weight = np.log2(total_samples + 1)
        final_score = proportion_diff * sample_weight
        
        return {
            'total_samples': total_samples,
            'meth_proportion': meth_proportion,
            'unmeth_proportion': unmeth_proportion,
            'proportion_diff': proportion_diff,
            'sample_weight': sample_weight,
            'final_score': final_score
        }
    
    def load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Load and process the methylation data from both haplotype files."""
        try:
            h1_df = pd.read_csv(self.h1_path, sep='\t')
            h2_df = pd.read_csv(self.h2_path, sep='\t')
            
            # Validate required columns
            required_cols = ['chrom', 'start', 'end', 'num_unmethylated', 'num_methylated', 
                           'unmethylated_samples', 'methylated_samples']
            for df in [h1_df, h2_df]:
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    raise ValueError(f"Missing required columns: {missing_cols}")
            
            # Add chromosome order for sorting
            chrom_order = {f'chr{i}': i for i in list(range(1, 23)) + ['X', 'Y']}
            
            for df in [h1_df, h2_df]:
                # Calculate methylation metrics
                metrics = df.apply(self.calculate_methylation_metrics, axis=1)
                
                # Add metrics as new columns
                df['total_samples'] = metrics.apply(lambda x: x['total_samples'])
                df['meth_proportion'] = metrics.apply(lambda x: x['meth_proportion'])
                df['unmeth_proportion'] = metrics.apply(lambda x: x['unmeth_proportion'])
                df['proportion_diff'] = metrics.apply(lambda x: x['proportion_diff'])
                df['sample_weight'] = metrics.apply(lambda x: x['sample_weight'])
                df['meth_score'] = metrics.apply(lambda x: x['final_score'])
                
                # Add score interpretation
                df['score_interpretation'] = df['meth_score'].apply(
                    lambda x: 'Strongly Methylated' if x > 1 else
                             'Moderately Methylated' if x > 0.5 else
                             'Weakly Methylated' if x > 0 else
                             'Neutral' if x == 0 else
                             'Weakly Unmethylated' if x > -0.5 else
                             'Moderately Unmethylated' if x > -1 else
                             'Strongly Unmethylated'
                )
                
                # Add chromosome order
                df['chrom_order'] = df['chrom'].map(lambda x: chrom_order.get(x, 999))
            
            return h1_df, h2_df
            
        except Exception as e:
            logger.error(f"Error loading data: {str(e)}")
            raise

    def format_metrics_for_tsv(self, df: pd.DataFrame) -> pd.DataFrame:
        """Format methylation metrics for TSV output."""
        df['score_calculation'] = df.apply(
            lambda row: (
                f"Total Samples: {row['total_samples']}\n"
                f"Methylated: {row['num_methylated']} ({row['meth_proportion']:.3f})\n"
                f"Unmethylated: {row['num_unmethylated']} ({row['unmeth_proportion']:.3f})\n"
                f"Proportion Difference: {row['proportion_diff']:.3f}\n"
                f"Sample Weight (log2(total + 1)): {row['sample_weight']:.3f}\n"
                f"Final Score: {row['proportion_diff']:.3f} × {row['sample_weight']:.3f} = {row['meth_score']:.3f}\n"
                f"Interpretation: {row['score_interpretation']}"
            ),
            axis=1
        )
        
        columns = [
            'chrom', 'start', 'end', 'meth_score', 'score_calculation',
            'num_methylated', 'num_unmethylated', 'total_samples',
            'meth_proportion', 'unmeth_proportion', 'proportion_diff',
            'sample_weight', 'score_interpretation',
            'methylated_samples', 'unmethylated_samples'
        ]
        
        return df[columns]

    def select_top_regions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Select regions with the most extreme methylation scores."""
        if len(df) < self.n_regions:
            logger.warning(f"Requested {self.n_regions} regions but only {len(df)} available")
            return df
            
        n_each = self.n_regions // 2
        top_methylated = df.nlargest(n_each, 'meth_score')
        top_unmethylated = df.nsmallest(n_each, 'meth_score')
        
        selected_regions = pd.concat([top_methylated, top_unmethylated])
        selected_regions = selected_regions.sort_values(['chrom_order', 'start']).reset_index(drop=True)
        
        return selected_regions
    
    def select_random_regions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Select random regions from the dataset."""
        if len(df) < self.n_regions:
            logger.warning(f"Requested {self.n_regions} regions but only {len(df)} available")
            return df.sample(n=len(df), random_state=42)
            
        return df.sample(n=self.n_regions, random_state=42).sort_values(['chrom_order', 'start']).reset_index(drop=True)
    
    def select_high_coverage_regions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Select regions with the highest total sample coverage."""
        if len(df) < self.n_regions:
            logger.warning(f"Requested {self.n_regions} regions but only {len(df)} available")
            return df
            
        high_coverage = df.nlargest(self.n_regions, 'total_samples')
        return high_coverage.sort_values(['chrom_order', 'start']).reset_index(drop=True)

    def create_static_plot(self, df: pd.DataFrame, plot_title: str, output_filename: str):
        """Create a static PNG plot for the selected regions."""
        # Create separate figures for each haplotype
        h1_subset = df[df['haplotype'] == 'H1']
        h2_subset = df[df['haplotype'] == 'H2']
        
        # Get max count for consistent scaling
        max_count = max(
            df['num_methylated'].max(),
            df['num_unmethylated'].max()
        )
        
        # Function to create a figure for one haplotype
        def create_haplotype_figure(data, title):
            # Calculate dynamic figure height
            fig_height = max(
                self.min_fig_height,
                len(data) * self.fig_height_per_region + 200  # Extra space for title and legend
            )
            
            fig = go.Figure()
            
            # Add unmethylated bars (left side)
            fig.add_trace(go.Bar(
                y=data['chrom'] + ':' + data['start'].astype(str) + '-' + data['end'].astype(str),
                x=-data['num_unmethylated'],
                orientation='h',
                name='Unmethylated',
                marker=dict(color=self.colors['unmethylated']),
                text=data['num_unmethylated'],
                textposition='inside',
            ))
            
            # Add methylated bars (right side)
            fig.add_trace(go.Bar(
                y=data['chrom'] + ':' + data['start'].astype(str) + '-' + data['end'].astype(str),
                x=data['num_methylated'],
                orientation='h',
                name='Methylated',
                marker=dict(color=self.colors['methylated']),
                text=data['num_methylated'],
                textposition='inside',
            ))
            
            # Update layout for better static rendering
            fig.update_layout(
                title=dict(
                    text=f"{plot_title} - {title}",
                    x=0.5,
                    xanchor='center',
                    y=0.95,
                    yanchor='top',
                    font=dict(size=16)
                ),
                xaxis_title="Sample Count",
                yaxis_title="Genomic Region",
                barmode='relative',
                bargap=0.2,
                bargroupgap=0.1,
                showlegend=True,
                legend=dict(
                    x=0.01,
                    y=-0.15,
                    orientation='h',
                    xanchor='left',
                    yanchor='top'
                ),
                margin=dict(
                    l=200,   # Left margin for region labels
                    r=50,    # Right margin
                    t=100,   # Top margin for title
                    b=100    # Bottom margin for legend
                ),
                plot_bgcolor='white',
                width=self.fig_width,
                height=fig_height,
                font=dict(size=12)  # Increase base font size for better readability
            )
            
            # Customize axes for static display
            fig.update_xaxes(
                range=[-max_count-0.5, max_count+0.5],
                zeroline=True,
                zerolinewidth=1,
                zerolinecolor='black',
                gridcolor='lightgrey',
                gridwidth=1,
                tickmode='array',
                tickvals=list(range(-max_count, max_count+1)),
                ticktext=[str(abs(x)) for x in range(-max_count, max_count+1)],
                tickangle=0
            )
            
            fig.update_yaxes(
                gridcolor='lightgrey',
                gridwidth=1,
                autorange="reversed"
            )
            
            # Add center line
            fig.add_vline(x=0, line_width=1, line_color="black")
            
            # Add score calculation explanation
            fig.add_annotation(
                text="Score = (methylated - unmethylated) / total × log2(total + 1)",
                xref="paper",
                yref="paper",
                x=0.99,
                y=-0.15,
                showarrow=False,
                font=dict(size=10, color="gray"),
                align="right"
            )
            
            return fig
        
        # Create and save plots for each haplotype
        base_filename = output_filename.replace('.png', '')
        
        # H1 plot
        fig_h1 = create_haplotype_figure(h1_subset, "Haplotype 1")
        fig_h1.write_image(str(self.output_dir / f"{base_filename}_H1.png"))
        
        # H2 plot
        fig_h2 = create_haplotype_figure(h2_subset, "Haplotype 2")
        fig_h2.write_image(str(self.output_dir / f"{base_filename}_H2.png"))
    
    def save_tsv(self, df: pd.DataFrame, filename: str):
        """Save the DataFrame to a TSV file."""
        df.to_csv(self.output_dir / filename, sep='\t', index=False)
    
    def run(self):
        """Execute the visualization pipeline."""
        try:
            logger.info("Loading data...")
            h1_df, h2_df = self.load_data()
            
            # Select and format regions for both haplotypes
            h1_top_regions = self.format_metrics_for_tsv(self.select_top_regions(h1_df))
            h2_top_regions = self.format_metrics_for_tsv(self.select_top_regions(h2_df))
            
            h1_random_regions = self.format_metrics_for_tsv(self.select_random_regions(h1_df))
            h2_random_regions = self.format_metrics_for_tsv(self.select_random_regions(h2_df))
            
            h1_high_coverage = self.format_metrics_for_tsv(self.select_high_coverage_regions(h1_df))
            h2_high_coverage = self.format_metrics_for_tsv(self.select_high_coverage_regions(h2_df))
            
            # Combine data for all region types
            top_regions_df = pd.concat([
                h1_top_regions.assign(haplotype='H1'),
                h2_top_regions.assign(haplotype='H2')
            ])
            
            random_regions_df = pd.concat([
                h1_random_regions.assign(haplotype='H1'),
                h2_random_regions.assign(haplotype='H2')
            ])
            
            high_coverage_df = pd.concat([
                h1_high_coverage.assign(haplotype='H1'),
                h2_high_coverage.assign(haplotype='H2')
            ])
            
            logger.info("Creating visualizations...")
            self.create_static_plot(top_regions_df, 'Top Discrepant Regions', 'top_discrepant_regions.png')
            self.create_static_plot(random_regions_df, 'Random Regions', 'random_regions.png')
            self.create_static_plot(high_coverage_df, 'Highest Coverage Regions', 'high_coverage_regions.png')
            
            logger.info("Saving TSV files...")
            self.save_tsv(top_regions_df, 'top_discrepant_regions.tsv')
            self.save_tsv(random_regions_df, 'random_regions.tsv')
            self.save_tsv(high_coverage_df, 'high_coverage_regions.tsv')
            
            logger.info("Visualization complete!")
            
        except Exception as e:
            logger.error(f"Error during visualization: {str(e)}")
            raise

def main():
    """Parse command line arguments and run the visualization."""
    parser = argparse.ArgumentParser(description="Visualize methylation patterns across haplotypes")
    parser.add_argument("--h1", required=True, help="Path to H1 haplotype BED file")
    parser.add_argument("--h2", required=True, help="Path to H2 haplotype BED file")
    parser.add_argument("--output", required=True, help="Output directory for plots")
    parser.add_argument("--n-regions", type=int, default=50, 
                       help="Number of regions to visualize for each category")
    
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