#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import matplotlib.gridspec as gridspec

# Genomic Features Constants
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
            # Read the file with headers
            df = pd.read_csv(bed_file, sep='\t')
            
            # Ensure required columns exist
            required_cols = ['chrom', 'start', 'end']
            if not all(col in df.columns for col in required_cols):
                print(f"Warning: Missing required columns in {bed_file.name}")
                continue
            
            # Convert start and end to numeric, replacing errors with NaN
            df['start'] = pd.to_numeric(df['start'], errors='coerce')
            df['end'] = pd.to_numeric(df['end'], errors='coerce')
            
            # Drop any rows where conversion failed
            df = df.dropna(subset=['start', 'end'])
            
            # Convert to integers
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            
            sample_name = bed_file.stem.split('_')[0]
            df['sample'] = sample_name
            
            # Use existing size column if present, otherwise calculate it
            if 'size' not in df.columns:
                df['size'] = df['end'] - df['start']
            
            # Only keep rows with positive size
            df = df[df['size'] > 0]
            
            if len(df) > 0:
                dfs.append(df)
            else:
                print(f"Warning: No valid regions found in {bed_file.name}")
                
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

def calculate_enhanced_bin_stats(data_dict, chrom, n_bins):
    """Calculate detailed bin statistics including counts and size variation"""
    bin_edges = np.linspace(0, CHROM_SIZES[chrom], n_bins + 1)
    bin_stats = []
    
    # Initialize bin statistics
    for bin_idx in range(n_bins):
        bin_stats.append({
            'bin': bin_idx,
            'position_mb': int(bin_edges[bin_idx] / 1e6),
            'counts_per_sample': [],
            'sizes': [],
            'mean_count': 0,
            'std_count': 0,
            'mean_size': 0,
            'std_size': 0,
            'median_size': 0,
            'size_range': (0, 0)
        })
    
    # Collect statistics for each sample
    for sample in sorted(data_dict.keys()):
        data = data_dict[sample]
        chrom_data = data[data['chrom'] == chrom]
        positions = (chrom_data['start'] + chrom_data['end']) / 2
        sizes = chrom_data['size'].values
        
        # Get bin indices for each position
        bin_indices = np.digitize(positions, bin_edges) - 1
        
        # Count regions and collect sizes for each bin
        for bin_idx in range(n_bins):
            mask = bin_indices == bin_idx
            count = np.sum(mask)
            bin_stats[bin_idx]['counts_per_sample'].append(count)
            if count > 0:
                bin_stats[bin_idx]['sizes'].extend(sizes[mask])
    
    # Calculate summary statistics for each bin
    for stats in bin_stats:
        if stats['counts_per_sample']:
            stats['mean_count'] = np.mean(stats['counts_per_sample'])
            stats['std_count'] = np.std(stats['counts_per_sample'])
        
        if stats['sizes']:
            stats['mean_size'] = np.mean(stats['sizes'])
            stats['std_size'] = np.std(stats['sizes'])
            stats['median_size'] = np.median(stats['sizes'])
            stats['size_range'] = (np.min(stats['sizes']), np.max(stats['sizes']))
    
    return bin_stats

def add_special_regions(ax, chrom):
    """Add special genomic regions to the plot"""
    if chrom in SPECIAL_REGIONS:
        for start, end, desc, region_type in SPECIAL_REGIONS[chrom]:
            color = REGION_COLORS[region_type]
            # Add shaded region
            ax.axvspan(start/1e6, end/1e6, color=color, alpha=0.2, 
                      label=f'{desc} ({region_type})')
            # Add boundary lines
            ax.axvline(x=start/1e6, color='gray', linestyle=':', linewidth=0.5)
            ax.axvline(x=end/1e6, color='gray', linestyle=':', linewidth=0.5)

def plot_enhanced_bin_stats(chrom, m_data_dict, u_data_dict, output_dir, condition='M'):
    """Create visualization showing both count and size statistics per bin"""
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10
    })
    
    # Create figure with two panels
    fig = plt.figure(figsize=(20, 12))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    fig.suptitle(f'{chrom} Region Statistics - {"Methylated" if condition=="M" else "Unmethylated"}',
                fontsize=16, y=0.95)
    
    data_dict = m_data_dict if condition == 'M' else u_data_dict
    
    # Calculate bin statistics with smaller bins for better resolution
    bin_size = 5_000_000  # 5Mb bins instead of 10Mb
    n_bins = int(np.ceil(CHROM_SIZES[chrom] / bin_size))
    bin_stats = calculate_enhanced_bin_stats(data_dict, chrom, n_bins)
    
    # Extract x-axis positions (bin centers in Mb)
    x_pos = [stats['position_mb'] for stats in bin_stats]
    
    # Colors and style
    fill_color = '#FF9999' if condition == 'M' else '#9999FF'
    line_color = '#CC0000' if condition == 'M' else '#0000CC'
    edge_color = '#800000' if condition == 'M' else '#000080'
    
    # Panel 1: Region Counts
    mean_counts = [stats['mean_count'] for stats in bin_stats]
    std_counts = [stats['std_count'] for stats in bin_stats]
    
    # Plot mean counts with improved styling
    bars1 = ax1.bar(x_pos, mean_counts, width=4, alpha=0.4, color=fill_color, 
                   edgecolor=edge_color, linewidth=1)
    ax1.errorbar(x_pos, mean_counts, yerr=std_counts, fmt='none', 
                color=line_color, capsize=3, capthick=1, elinewidth=1)
    
    # Add count labels more selectively
    max_count = max(mean_counts)
    threshold = max_count * 0.1  # Only label bars that are at least 10% of max height
    for i, (mean, std) in enumerate(zip(mean_counts, std_counts)):
        if mean > threshold:
            ax1.text(x_pos[i], mean + std + (max_count * 0.02), 
                    f'{mean:.0f}±{std:.0f}',
                    ha='center', va='bottom', rotation=45, fontsize=8)
    
    ax1.set_xlabel('Position (Mb)')
    ax1.set_ylabel('Average Number of Regions per Sample')
    ax1.set_title('Region Count Distribution', pad=20)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Panel 2: Size Statistics
    mean_sizes = [stats['mean_size'] for stats in bin_stats]
    std_sizes = [stats['std_size'] for stats in bin_stats]
    median_sizes = [stats['median_size'] for stats in bin_stats]
    
    # Plot size statistics with violin or box plots
    size_data = [stats['sizes'] for stats in bin_stats]
    violin_parts = ax2.violinplot(size_data, positions=x_pos, 
                                showmeans=True, showextrema=True)
    
    # Customize violin plot colors
    for pc in violin_parts['bodies']:
        pc.set_facecolor(fill_color)
        pc.set_alpha(0.3)
    violin_parts['cmeans'].set_color(edge_color)
    violin_parts['cmins'].set_color(edge_color)
    violin_parts['cmaxes'].set_color(edge_color)
    violin_parts['cbars'].set_color(edge_color)
    
    # Add size labels more selectively
    max_size = max(mean_sizes)
    threshold = max_size * 0.1
    for i, (mean, median) in enumerate(zip(mean_sizes, median_sizes)):
        if mean > threshold:
            ax2.text(x_pos[i], mean + (max_size * 0.02),
                    f'Mean: {int(mean):,}\nMedian: {int(median):,}',
                    ha='center', va='bottom', rotation=45, fontsize=8)
    
    ax2.set_xlabel('Position (Mb)')
    ax2.set_ylabel('Region Size (bp)')
    ax2.set_title('Region Size Distribution', pad=20)
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    # Add genomic features to both panels
    for ax in [ax1, ax2]:
        add_special_regions(ax, chrom)
        
        if chrom in CENTROMERES:
            cent_start, cent_end = CENTROMERES[chrom]
            ax.axvspan(cent_start/1e6, cent_end/1e6,
                      color='gray', alpha=0.2,
                      label='Centromere')
        
        # Add legend
        ax.legend(loc='upper right', bbox_to_anchor=(1, 1),
                 ncol=1, fontsize=8)
        
        # Set reasonable y-axis limits
        if ax == ax1:
            ax.set_ylim(0, max_count * 1.2)
        else:
            ax.set_ylim(0, max_size * 1.2)
    
    plt.tight_layout()
    output_file = output_dir / f'{chrom}_{condition}_bin_statistics.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return bin_stats

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate detailed region size statistics and visualizations')
    parser.add_argument('--base_dir', required=True, help='Base directory containing condition folders')
    parser.add_argument('--output_dir', required=True, help='Output directory for plots and statistics')
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Starting region size analysis...")
    
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
    
    # Create summary directory for combined statistics
    summary_dir = output_dir / 'summary'
    summary_dir.mkdir(exist_ok=True)
    
    # Generate global summary statistics
    print("\nGenerating global statistics...")
    global_stats = []
    for condition, data_dict in [('Methylated', m_data), ('Unmethylated', u_data)]:
        for sample, data in data_dict.items():
            stats = {
                'condition': condition,
                'sample': sample,
                'total_regions': len(data),
                'mean_size': data['size'].mean(),
                'median_size': data['size'].median(),
                'std_size': data['size'].std(),
                'min_size': data['size'].min(),
                'max_size': data['size'].max()
            }
            # Add chromosome-specific counts
            for chrom in data['chrom'].unique():
                stats[f'{chrom}_count'] = len(data[data['chrom'] == chrom])
            global_stats.append(stats)
    
    # Save global statistics
    global_stats_df = pd.DataFrame(global_stats)
    global_stats_file = summary_dir / 'global_statistics.csv'
    global_stats_df.to_csv(global_stats_file, index=False)
    print(f"Saved global statistics to {global_stats_file}")
    
    # Process each chromosome
    chromosomes = sorted([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'],
                        key=lambda x: int(x.replace('chr', '')) 
                        if x.replace('chr', '').isdigit() else float('inf'))
    
    # Create chromosome-specific directories
    chrom_stats = []
    for chrom in chromosomes:
        print(f"\nProcessing {chrom}...")
        chrom_dir = output_dir / chrom
        chrom_dir.mkdir(exist_ok=True)
        
        # Generate methylated and unmethylated plots
        m_stats = plot_enhanced_bin_stats(chrom, m_data, u_data, chrom_dir, 'M')
        u_stats = plot_enhanced_bin_stats(chrom, m_data, u_data, chrom_dir, 'U')
        
        # Collect chromosome-level statistics
        for condition, stats in [('M', m_stats), ('U', u_stats)]:
            for bin_stat in stats:
                if bin_stat['mean_count'] > 0:  # Only include bins with data
                    chrom_stats.append({
                        'chromosome': chrom,
                        'condition': 'Methylated' if condition == 'M' else 'Unmethylated',
                        'bin_position_mb': bin_stat['position_mb'],
                        'mean_regions_per_sample': bin_stat['mean_count'],
                        'std_regions': bin_stat['std_count'],
                        'mean_region_size': bin_stat['mean_size'],
                        'median_region_size': bin_stat['median_size'],
                        'min_region_size': bin_stat['size_range'][0],
                        'max_region_size': bin_stat['size_range'][1]
                    })
    
    # Save chromosome-level statistics
    chrom_stats_df = pd.DataFrame(chrom_stats)
    chrom_stats_file = summary_dir / 'chromosome_bin_statistics.csv'
    chrom_stats_df.to_csv(chrom_stats_file, index=False)
    print(f"Saved chromosome statistics to {chrom_stats_file}")
    
    print("\nAnalysis complete!")
    print("\nSummary of outputs:")
    print(f"1. Global statistics: {global_stats_file}")
    print(f"2. Chromosome bin statistics: {chrom_stats_file}")
    print(f"3. Individual chromosome plots and data: {output_dir}")
    
    # Print some key findings
    print("\nKey Statistics:")
    for condition in ['Methylated', 'Unmethylated']:
        condition_stats = global_stats_df[global_stats_df['condition'] == condition]
        print(f"\n{condition} Regions:")
        print(f"  Average regions per sample: {condition_stats['total_regions'].mean():.1f}")
        print(f"  Mean region size: {condition_stats['mean_size'].mean():.1f} bp")
        print(f"  Size range: {condition_stats['min_size'].min():.0f} - {condition_stats['max_size'].max():.0f} bp")

if __name__ == "__main__":
    main()