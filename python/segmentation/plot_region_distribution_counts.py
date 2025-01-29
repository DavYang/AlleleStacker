#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import matplotlib.gridspec as gridspec

# Updated centromere positions with ranges (hg38)
CENTROMERES = {
    'chr1': (122026459, 124932724),
    'chr2': (92188145, 94090557),
    'chr3': (90772458, 93655574),
    'chr4': (49712061, 51743951),
    'chr5': (46485900, 50059807),
    'chr6': (58553888, 59829934),
    'chr7': (58169653, 61528020),
    'chr8': (44033744, 45877265),
    'chr9': (43389635, 45518558),
    'chr10': (39686682, 41593521),
    'chr11': (51078348, 54425074),
    'chr12': (34769407, 37185252),
    'chr13': (16000000, 18051248),
    'chr14': (16000000, 18173523),
    'chr15': (17083673, 19725254),
    'chr16': (36311158, 38265669),
    'chr17': (22813679, 26616164),
    'chr18': (15460899, 20813083),
    'chr19': (24498980, 27190874),
    'chr20': (26436232, 30038348),
    'chr21': (10864560, 12915808),
    'chr22': (12954788, 15054318),
    'chrX': (58605579, 62412542),
    'chrY': (10316944, 10544039)
}

# Chromosome sizes (hg38)
CHROM_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
    'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

# Add the special genomic regions
SPECIAL_REGIONS = {
    # Acrocentric Chromosome Short Arms
    "chr13": [(0, 16000000, "Short arm (p)", "Acrocentric")],
    "chr14": [(0, 16000000, "Short arm (p)", "Acrocentric")],
    "chr15": [(0, 17000000, "Short arm (p)", "Acrocentric")],
    "chr21": [(0, 12000000, "Short arm (p)", "Acrocentric")],
    "chr22": [(0, 14000000, "Short arm (p)", "Acrocentric")],
    # Telomeric Regions and PAR
    "chrX": [
        (10000, 2781479, "PAR1", "Pseudoautosomal"),
        (155701383, 156030895, "PAR2", "Pseudoautosomal")
    ],
    "chrY": [
        (10000, 2781479, "PAR1", "Pseudoautosomal"),
        (56887903, 57217415, "PAR2", "Pseudoautosomal"),
        (13000000, 59373566, "Heterochromatin", "Heterochromatic"),
        (14000000, 15000000, "AZFa", "Ampliconic"),
        (18000000, 24000000, "AZFb", "Ampliconic"),
        (24000000, 28000000, "AZFc", "Ampliconic")
    ],
    # Major Heterochromatic Regions
    "chr1": [(121700000, 125100000, "Heterochromatin", "Heterochromatic")],
    "chr9": [(42200000, 45500000, "Heterochromatin", "Heterochromatic")],
    "chr16": [(35335801, 36185801, "Heterochromatin", "Heterochromatic")]
}

# Define colors for different region types
REGION_COLORS = {
    "Acrocentric": "#FFB6C1",  # Light pink
    "Pseudoautosomal": "#98FB98",  # Pale green
    "Heterochromatic": "#DDA0DD",  # Plum
    "Ampliconic": "#87CEEB"  # Sky blue
}

def load_bed_files(base_dir, condition):
    """Load and combine BED files for a given condition"""
    print(f"\nProcessing condition: {condition}")
    path = Path(base_dir) / condition
    dfs = []
    
    if not path.exists():
        print(f"Warning: Path {path} does not exist")
        return None
        
    bed_files = list(path.glob('*.bed'))
    print(f"Found {len(bed_files)} BED files")
    
    for bed_file in bed_files:
        print(f"Loading {bed_file.name}...")
        try:
            df = pd.read_csv(bed_file, sep='\t')
            sample_name = bed_file.stem.split('_')[0]
            df['sample'] = sample_name
            dfs.append(df)
        except Exception as e:
            print(f"Error loading {bed_file}: {str(e)}")
            continue
    
    if not dfs:
        print("No valid BED files found")
        return None
        
    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Total regions loaded: {len(combined_df)}")
    return combined_df

def filter_chromosomes(df):
    """Filter for standard chromosomes only"""
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    filtered_df = df[df['chrom'].isin(valid_chroms)]
    print(f"Filtered from {len(df)} to {len(filtered_df)} regions")
    return filtered_df

def calculate_density(data_dict, chrom, n_bins):
    """Calculate region density for all samples"""
    bin_edges = np.linspace(0, CHROM_SIZES[chrom], n_bins + 1)
    density_matrix = []
    sample_names = []
    
    # Sort samples to ensure consistent ordering
    for sample in sorted(data_dict.keys()):
        data = data_dict[sample]
        chrom_data = data[data['chrom'] == chrom]
        positions = (chrom_data['start'] + chrom_data['end']) / 2
        hist, _ = np.histogram(positions, bins=bin_edges)
        density_matrix.append(hist)
        sample_names.append(sample)
    
    return np.array(density_matrix), sample_names

def add_special_regions(ax, chrom):
    """Add special genomic regions to the plot"""
    if chrom in SPECIAL_REGIONS:
        for start, end, desc, region_type in SPECIAL_REGIONS[chrom]:
            color = REGION_COLORS[region_type]
            # Add shaded region
            ax.axvspan(start, end, color=color, alpha=0.2, 
                      label=f'{desc} ({region_type})')
            # Add boundary lines
            ax.axvline(x=start, color='gray', linestyle=':', linewidth=0.5)
            ax.axvline(x=end, color='gray', linestyle=':', linewidth=0.5)

def plot_chromosome_panels(chrom, m_data_dict, u_data_dict, output_dir, condition='M'):
    """Create three-panel visualization for a chromosome"""
    print(f"\nProcessing {condition} regions for chromosome {chrom}")
    
    # Increase font sizes
    plt.rcParams.update({'font.size': 12})  # Base font size increased by 15%
    
    # Create figure with three panels - heatmap first
    fig = plt.figure(figsize=(20, 15))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1.2, 1, 2], hspace=0.4)
    
    # Get data dictionary based on condition
    data_dict = m_data_dict if condition == 'M' else u_data_dict
    
    # Calculate bin edges and centers
    # bin_size = 10_000_000  # 10Mb bins 
    # bin_size = 1_000_000  # 1Mb bins
    bin_size = 100_000  # 100kb bins
    n_bins = int(np.ceil(CHROM_SIZES[chrom] / bin_size))
    bin_edges = np.arange(0, CHROM_SIZES[chrom] + bin_size, bin_size)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Panel 1: Heatmap (now first)
    ax1 = fig.add_subplot(gs[0])
    density_matrix, sample_names = calculate_density(data_dict, chrom, n_bins)
    
    # Use raw counts for visualization
    vmax = np.percentile(density_matrix, 95)  # Cap at 95th percentile for better contrast
    cmap = 'YlOrRd' if condition == 'M' else 'YlGnBu'
    
    im = ax1.imshow(density_matrix, aspect='auto', 
                    extent=[0, CHROM_SIZES[chrom], len(sample_names), 0],
                    cmap=cmap, vmin=0, vmax=vmax)
    
    ax1.set_yticks(np.arange(len(sample_names)) + 0.5)
    ax1.set_yticklabels(sample_names, fontsize=11)
    ax1.set_xlabel('Position (Mb)', fontsize=12)
    ax1.set_ylabel('Samples', fontsize=12)
    
    cbar = plt.colorbar(im, ax=ax1)
    cbar.set_label('Number of Regions', fontsize=12)  # Changed label to reflect raw counts
    
    # Panel 2: Consensus View
    ax2 = fig.add_subplot(gs[1])
    mean_density = np.mean(density_matrix, axis=0)
    std_density = np.std(density_matrix, axis=0)
    
    fill_color = 'red' if condition == 'M' else 'blue'
    line_color = 'darkred' if condition == 'M' else 'darkblue'
    
    ax2.fill_between(bin_centers, mean_density - std_density, 
                     mean_density + std_density, alpha=0.3,
                     color=fill_color)
    ax2.plot(bin_centers, mean_density, 
             color=line_color,
             label='Mean distribution')
    
    ax2.set_xlabel('Position (Mb)', fontsize=12)
    ax2.set_ylabel('Region Count', fontsize=12)
    
    # Panel 3: Stacked Distribution
    ax3 = fig.add_subplot(gs[2])
    
    all_samples = sorted(data_dict.keys())
    positions_list = []
    labels = []
    
    for sample in all_samples:
        data = data_dict[sample]
        chrom_data = data[data['chrom'] == chrom]
        if len(chrom_data) > 0:
            positions = (chrom_data['start'] + chrom_data['end']) / 2
            mean_size = int(chrom_data['size'].mean())
            positions_list.append(positions)
            labels.append(f'{sample} (n={len(chrom_data):,}, mean={mean_size:,}bp)')
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(positions_list)))
    ax3.hist(positions_list, bins=n_bins, stacked=True, 
             label=labels, color=colors, alpha=0.6)
    
    ax3.set_xlabel('Position (Mb)', fontsize=12)
    ax3.set_ylabel('Cumulative Region Count', fontsize=12)
    
    # Add special regions and format all panels
    for ax in [ax1, ax2, ax3]:
        add_special_regions(ax, chrom)
        
        if chrom in CENTROMERES:
            cent_start, cent_end = CENTROMERES[chrom]
            ax.axvspan(cent_start, cent_end, 
                      color='gray', alpha=0.2, 
                      label='Centromere' if ax == ax1 else "")
        
        ax.set_xlim(0, CHROM_SIZES[chrom])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda x, p: f'{int(x/1e6)}'))
    
    # Customize titles with larger font
    ax1.set_title('Region Density Heatmap by Sample', fontsize=14, pad=20)
    ax2.set_title(f'{"Methylated" if condition=="M" else "Unmethylated"} Consensus Distribution - {chrom}', 
                  fontsize=14, pad=20)
    ax3.set_title('Sample Distribution Breakdown', fontsize=14, pad=20)
    
    # Add legends with increased size
    # For special regions (first panel)
    handles, labels = ax1.get_legend_handles_labels()
    if handles:
        ax1.legend(handles, labels, bbox_to_anchor=(1.15, 1), 
                  loc='upper left', fontsize=11)
    
    # For sample distribution (third panel)
    ax3.legend(bbox_to_anchor=(1.15, 1), loc='upper left', 
              fontsize=11, title='Sample Information')
    
    plt.tight_layout()
    output_file = output_dir / f'{chrom}_{condition}_distribution.png'
    print(f"Saving plot to {output_file}")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate methylation distribution plots')
    parser.add_argument('--base_dir', required=True, help='Base directory containing condition folders')
    parser.add_argument('--output_dir', required=True, help='Output directory for plots')
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Starting visualization...")
    
    # Load all data
    conditions = ['H1_M', 'H2_M', 'H1_U', 'H2_U']
    all_data = {}
    
    for condition in conditions:
        df = load_bed_files(args.base_dir, condition)
        if df is not None:
            df = filter_chromosomes(df)
            all_data[condition] = df
    
    # Organize data by methylation status
    m_data = {}  # methylated data by sample
    u_data = {}  # unmethylated data by sample
    
    for condition, df in all_data.items():
        for sample in df['sample'].unique():
            sample_data = df[df['sample'] == sample]
            if condition.endswith('_M'):
                m_data[sample] = sample_data
            else:
                u_data[sample] = sample_data
    
    # Process each chromosome
    chromosomes = sorted([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'],
                        key=lambda x: int(x.replace('chr', '')) 
                        if x.replace('chr', '').isdigit() else float('inf'))
    
    for chrom in chromosomes:
        # Generate methylated and unmethylated plots
        plot_chromosome_panels(chrom, m_data, u_data, output_dir, 'M')
        plot_chromosome_panels(chrom, m_data, u_data, output_dir, 'U')
    
    print("\nVisualization complete!")

if __name__ == "__main__":
    main()