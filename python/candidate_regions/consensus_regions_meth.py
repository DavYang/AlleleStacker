#!/usr/bin/env python
import os
import argparse
from collections import defaultdict, namedtuple
import numpy as np
from typing import List, Dict, Tuple, Optional, Iterator
import csv
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from functools import partial

# Define the Event namedtuple
Event = namedtuple('Event', ['position', 'is_start', 'region'])

CHROM_ORDER = {f'chr{i}': i for i in range(1, 23)}
CHROM_ORDER.update({'chrX': 23, 'chrY': 24, 'chrM': 25})

def log_time(message: str):
    """Print message with timestamp"""
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{current_time}] {message}")

class Region:
    def __init__(self, chrom: str, start: int, end: int, sample_id: str, 
                 summary_label: Optional[str] = None, size: Optional[int] = None):
        if not isinstance(chrom, str) or not chrom.strip():
            raise ValueError(f"Chromosome must be a non-empty string, got {repr(chrom)}")
        if not isinstance(start, int) or start < 0:
            raise ValueError(f"Start position must be a non-negative integer, got {repr(start)}")
        if not isinstance(end, int) or end < 0:
            raise ValueError(f"End position must be a non-negative integer, got {repr(end)}")
        if end <= start:
            raise ValueError(f"End position ({end}) must be greater than start position ({start})")
            
        self.chrom = chrom
        self.start = start
        self.end = end
        self.sample_id = {sample_id} if isinstance(sample_id, str) else set(sample_id)
        self.summary_label = {summary_label} if summary_label else set()
        self.size = size or end - start
        self.contributing_samples = self.sample_id.copy()  # Track all contributing samples

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return CHROM_ORDER.get(self.chrom, 999) < CHROM_ORDER.get(other.chrom, 999)
        return (self.start, self.end) < (other.start, other.end)

    def __repr__(self):
        return f"Region(chrom='{self.chrom}', start={self.start}, end={self.end}, samples={self.sample_id})"

def calculate_size_stats(regions: List[Region]) -> Tuple[float, float]:
    """Calculate mean and standard deviation of region sizes"""
    if not regions:
        return 0.0, 0.0
    sizes = np.array([r.size for r in regions], dtype=np.float64)
    return np.mean(sizes), np.std(sizes)

def region_generator(file_path: str, sample_id: str) -> Iterator[Region]:
    """Generate Region objects from input file"""
    try:
        log_time(f"Reading regions from {file_path}")
        region_count = 0
        
        # Use a larger buffer size for faster file reading
        with open(file_path, 'r', buffering=1024*1024) as f:
            next(f)  # Skip header
            # Process lines in chunks for better performance
            for chunk in iter(lambda: f.readlines(1024*1024), []):
                for line in chunk:
                    try:
                        # Split only the needed fields for better performance
                        fields = line.strip().split('\t', 5)
                        chrom, start, end, summary_label, size = fields[:5]
                        yield Region(chrom, int(start), int(end), sample_id, summary_label, int(size))
                        region_count += 1
                        if region_count % 50000 == 0:  # Reduced logging frequency
                            log_time(f"Processed {region_count} regions from {file_path}")
                    except ValueError as e:
                        continue
                    
        log_time(f"Completed reading {region_count} regions from {file_path}")
    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")
        return []

def merge_regions(regions: List[Region], max_gap: int = 500, min_samples: int = 10) -> List[Region]:
    """Optimized merge_regions function"""
    if not regions:
        return []
    
    log_time(f"Starting region merging for {len(regions)} regions")
    
    # Pre-sort regions by chromosome and start position
    regions.sort(key=lambda r: (CHROM_ORDER.get(r.chrom, 999), r.start, r.end))
    
    # Group regions by chromosome using dict comprehension
    chrom_regions = {}
    current_chrom = regions[0].chrom
    start_idx = 0
    
    merged_regions = []
    
    for i, region in enumerate(regions):
        if region.chrom != current_chrom:
            chrom_regions[current_chrom] = regions[start_idx:i]
            current_chrom = region.chrom
            start_idx = i
    chrom_regions[current_chrom] = regions[start_idx:]
    
    for chrom, chr_regions in chrom_regions.items():
        log_time(f"Processing chromosome {chrom}")
        
        # Use sets for faster lookups
        sample_contributions = defaultdict(set)
        
        i = 0
        while i < len(chr_regions):
            current_region = chr_regions[i]
            
            # Check sample availability using set operations
            samples_available = all(
                not any(
                    prev_start <= current_region.end + max_gap and 
                    prev_end >= current_region.start - max_gap
                    for prev_start, prev_end in sample_contributions[sample]
                )
                for sample in current_region.sample_id
            )
            
            if not samples_available:
                i += 1
                continue
            
            # Rest of the merging logic remains similar but uses set operations
            current_end = current_region.end
            overlapping_regions = [current_region]
            contributing_samples = {sample: {current_region} 
                                 for sample in current_region.sample_id}
            
            j = i + 1
            while j < len(chr_regions):
                next_region = chr_regions[j]
                
                if next_region.start > current_end + max_gap:
                    break
                
                # Use set operations for sample validation
                new_samples = next_region.sample_id - contributing_samples.keys()
                sample_ok = len(new_samples) == len(next_region.sample_id)
                
                if sample_ok:
                    overlapping_regions.append(next_region)
                    current_end = max(current_end, next_region.end)
                    for sample in next_region.sample_id:
                        contributing_samples[sample] = {next_region}
                
                j += 1
            
            if len(contributing_samples) >= min_samples:
                consensus_region = Region(
                    chrom,
                    min(r.start for r in overlapping_regions),
                    max(r.end for r in overlapping_regions),
                    next(iter(contributing_samples.keys())),
                    None
                )
                consensus_region.sample_id = set(contributing_samples.keys())
                consensus_region.contributing_samples = set().union(*(r.sample_id for r in overlapping_regions))
                merged_regions.append(consensus_region)
                
                # Update sample contributions using set operations
                for sample, regions in contributing_samples.items():
                    for region in regions:
                        sample_contributions[sample].add((region.start, region.end))
            
            i = j
            
        log_time(f"Created {len(merged_regions)} consensus regions for chromosome {chrom}")
    
    log_time(f"Merging completed: {len(merged_regions)} total consensus regions")
    return merged_regions

def write_consensus_regions(consensus_regions: List[Region], output_file: str) -> None:
    """Write consensus regions to output file with improved sample tracking"""
    log_time(f"Writing {len(consensus_regions)} consensus regions to {output_file}")
    
    with open(output_file, 'w', buffering=8192) as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            "chrom", "start", "end", "num_samples", "samples", 
            "summary_labels", "size", "contributing_samples"
        ])
        
        for idx, region in enumerate(sorted(consensus_regions)):
            writer.writerow([
                region.chrom,
                region.start,
                region.end,
                len(region.sample_id),
                ",".join(sorted(region.sample_id)),
                ",".join(sorted(region.summary_label)),
                region.size,
                ",".join(sorted(region.contributing_samples))
            ])
            
            if (idx + 1) % 10000 == 0:
                log_time(f"Wrote {idx + 1} regions")
    
    log_time(f"Completed writing consensus regions")

def process_sample(sample: str, input_dir: str, haplotype: str) -> Tuple[str, List[Region]]:
    """Process a single sample"""
    log_time(f"Processing sample {sample} for haplotype {haplotype}")
    
    sample_file = os.path.join(input_dir, f"{sample}_{haplotype}_M.bed")
    if not os.path.exists(sample_file):
        log_time(f"Warning: File not found for sample {sample}")
        return sample, []
    
    regions = list(region_generator(sample_file, sample))
    log_time(f"Completed processing {sample}: {len(regions)} regions")
    return sample, regions

def process_haplotype(input_dir: str, output_file: str, max_gap: int, 
                     min_samples: int, sample_list: List[str], haplotype: str) -> Dict:
    """Process regions for a specific haplotype with parallel processing"""
    log_time(f"Starting haplotype {haplotype} processing")
    
    sample_stats = {}
    all_regions = []
    
    # Increase max_workers for more parallelism
    max_workers = min(64, len(sample_list))  # Increased from 32
    log_time(f"Processing {len(sample_list)} samples using {max_workers} workers")
    
    # Use larger chunk sizes for better thread utilization
    chunk_size = max(1, len(sample_list) // (max_workers * 2))
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        process_fn = partial(process_sample, input_dir=input_dir, haplotype=haplotype)
        results = list(executor.map(process_fn, sample_list, chunksize=chunk_size))
        
        # Use list comprehension for faster processing
        all_regions = [region for _, regions in results for region in regions if regions]
        sample_stats = {
            sample: {
                'count': len(regions),
                **dict(zip(['mean_size', 'size_std'], calculate_size_stats(regions)))
            }
            for sample, regions in results if regions
        }

    if not all_regions:
        log_time(f"No valid regions found in {input_dir}")
        return sample_stats

    log_time(f"Processing {len(all_regions)} total regions")
    consensus_regions = merge_regions(all_regions, max_gap=max_gap, min_samples=min_samples)
    
    # Validate consensus regions
    valid_consensus = []
    for region in consensus_regions:
        if len(region.contributing_samples) >= min_samples:
            valid_consensus.append(region)
        else:
            log_time(f"Warning: Removing consensus region with insufficient samples: {region}")
    
    write_consensus_regions(valid_consensus, output_file)
    
    mean_size, std_size = calculate_size_stats(valid_consensus)
    sample_stats['consensus'] = {
        'count': len(valid_consensus),
        'mean_size': mean_size,
        'size_std': std_size
    }
    
    log_time(f"Completed haplotype {haplotype} processing")
    return sample_stats

def main():
    log_time("Starting consensus region generation")
    
    parser = argparse.ArgumentParser(description="Generate consensus methylated regions for H1 and H2 haplotypes.")
    parser.add_argument("--input_dir", required=True, help="Input directory containing H1_M and H2_M subdirectories")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output BED files")
    parser.add_argument("--output_dir", required=True, help="Output directory for consensus BED files")
    parser.add_argument("--sample-list-file", required=True, help="File containing sample IDs")
    parser.add_argument("--max-gap", type=int, default=5000, help="Maximum gap between regions")
    parser.add_argument("--min-samples", type=int, default=2, help="Minimum samples required")
    args = parser.parse_args()

    log_time("Creating output directory")
    os.makedirs(args.output_dir, exist_ok=True)

    log_time("Reading sample list")
    with open(args.sample_list_file) as f:
        sample_list = [line.strip() for line in f if line.strip()]
    log_time(f"Found {len(sample_list)} samples")

    stats = {}
    for hap in ['H1', 'H2']:
        log_time(f"\nProcessing {hap} haplotype...")
        input_dir = os.path.join(args.input_dir, f"{hap}_M")
        output_file = os.path.join(args.output_dir, f"{args.output_prefix}_{hap}.bed")
        stats[hap] = process_haplotype(input_dir, output_file, args.max_gap, 
                                     args.min_samples, sample_list, hap)

    log_time("Writing statistics file")
    stats_file = os.path.join(args.output_dir, f"{args.output_prefix}_statistics.tsv")
    with open(stats_file, 'w', buffering=8192) as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            "Sample", "Haplotype", "Count", 
            "Mean_Size", "Size_StdDev"
        ])
        for hap in ['H1', 'H2']:
            for sample, stats_dict in sorted(stats[hap].items()):
                writer.writerow([
                    sample,
                    hap,
                    stats_dict['count'],
                    f"{stats_dict['mean_size']:.2f}",
                    f"{stats_dict['size_std']:.2f}"
                ])

    # Calculate and print summary statistics
    log_time("\nSummary Statistics:")
    for hap in ['H1', 'H2']:
        sample_sizes = np.array([
            stats[hap][sample]['mean_size'] 
            for sample in stats[hap] 
            if sample != 'consensus'
        ])
        sample_counts = np.array([
            stats[hap][sample]['count']
            for sample in stats[hap]
            if sample != 'consensus'
        ])
        consensus_stats = stats[hap]['consensus']
        
        print(f"\n{hap} Statistics:")
        print(f"Average number of regions per sample: {np.mean(sample_counts):.2f} ± {np.std(sample_counts):.2f}")
        print(f"Average region size across samples: {np.mean(sample_sizes):.2f} ± {np.std(sample_sizes):.2f}")
        print(f"Consensus regions: {consensus_stats['count']} regions")
        print(f"Consensus region size: {consensus_stats['mean_size']:.2f} ± {consensus_stats['size_std']:.2f}")

    log_time(f"\nProcessing complete. Output files written to {args.output_dir}")

if __name__ == "__main__":
    main()