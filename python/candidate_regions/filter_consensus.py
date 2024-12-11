#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict
import bisect
from typing import List, Dict, Set, Tuple
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import multiprocessing
from tqdm import tqdm

logging.basicConfig(
    level=logging.DEBUG,  # Change to DEBUG level
    format='%(asctime)s - %(levelname)s - %(message)s'
)

CHROM_ORDER = {f'chr{i}': i for i in range(1, 23)}
CHROM_ORDER.update({'chrX': 23, 'chrY': 24, 'chrM': 25})

class RegionIndex:
    def __init__(self):
        self.regions: Dict[str, List[Tuple[int, int]]] = {}
        self.sorted = False
        self.start_positions: Dict[str, List[int]] = {}
    
    def add_regions(self, chrom: str, regions: List[Tuple[int, int]]):
        if chrom not in self.regions:
            self.regions[chrom] = []
        self.regions[chrom].extend(regions)
        self.sorted = False
    
    def sort_all(self):
        for chrom in self.regions:
            self.regions[chrom].sort()
            self.start_positions[chrom] = [start for start, _ in self.regions[chrom]]
        self.sorted = True

    def find_containing_region(self, chrom: str, target_start: int, target_end: int) -> bool:
        """
        Check if any region in our index completely contains the target region.
        A region contains the target if: start <= target_start AND end >= target_end
        """
        if not self.sorted:
            self.sort_all()
            
        if chrom not in self.regions:
            return False
            
        # Find insertion point for target_start
        idx = bisect.bisect_right(self.start_positions[chrom], target_start)
        
        # Check region before insertion point
        if idx > 0:
            prev_start = self.start_positions[chrom][idx - 1]
            prev_end = self.regions[chrom][idx - 1][1]
            
            if prev_start <= target_start and prev_end >= target_end:
                return True
                
        # Check region at insertion point if it exists
        if idx < len(self.regions[chrom]):
            curr_start = self.start_positions[chrom][idx]
            curr_end = self.regions[chrom][idx][1]
            
            if curr_start <= target_start and curr_end >= target_end:
                return True
                
        return False
      
    def find_partial_overlap(self, chrom: str, target_start: int, target_end: int) -> bool:
        """
        Check if any region in our index partially overlaps the target region.
        """
        if not self.sorted:
            self.sort_all()
            
        if chrom not in self.regions:
            return False
            
        # Find insertion point for target_start
        idx = bisect.bisect_right(self.start_positions[chrom], target_start)
        
        # Check region before insertion point
        if idx > 0:
            prev_start = self.start_positions[chrom][idx - 1]
            prev_end = self.regions[chrom][idx - 1][1]
            
            # Check for partial overlap conditions
            if (prev_start <= target_start < prev_end or
                prev_start < target_end <= prev_end or
                target_start < prev_start < target_end):
                return True
                
        # Check region at insertion point if it exists
        if idx < len(self.regions[chrom]):
            curr_start = self.start_positions[chrom][idx]
            curr_end = self.regions[chrom][idx][1]
            
            # Check for partial overlap conditions
            if (curr_start <= target_start < curr_end or
                curr_start < target_end <= curr_end or
                target_start < curr_start < target_end):
                return True
                
        return False

def process_sample_file(args: Tuple[str, str]) -> Tuple[str, Dict[str, List[Tuple[int, int]]]]:
    """Process a single sample file and return its regions by chromosome."""
    sample_file, sample_id = args
    regions_by_chrom = defaultdict(list)
    
    try:
        with open(sample_file, 'r') as f:
            first_line = True
            for line in f:
                if first_line:
                    first_line = False
                    continue
                if line.startswith('#'): continue
                
                fields = line.strip().split()
                if len(fields) >= 3:
                    try:
                        chrom, start, end = fields[:3]
                        regions_by_chrom[chrom].append((int(start), int(end)))
                    except ValueError:
                        continue
    except Exception as e:
        logging.error(f"Error reading file {sample_file}: {str(e)}")
        return sample_id, {}
    
    return sample_id, dict(regions_by_chrom)

def read_consensus_regions(file_path: str) -> Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, Set[str]]]:
    """Read consensus regions from a BED file with sample information."""
    regions = defaultdict(list)
    region_samples = {}
    first_line = True
    
    logging.info(f"Reading consensus regions from {file_path}")
    with open(file_path, 'r') as f:
        for line in f:
            if first_line:
                first_line = False
                continue
            if line.startswith('#'): continue
            
            fields = line.strip().split('\t')
            if len(fields) >= 7:  # Must have chrom, start, end, num_samples, samples, labels, size
                chrom, start, end = fields[:3]
                samples = set(fields[4].split(',')) if fields[4] != "." else set()
                try:
                    region_key = (int(start), int(end))
                    regions[chrom].append(region_key)
                    region_samples[f"{chrom}:{start}-{end}"] = samples
                except ValueError:
                    continue
    
    # Sort regions by chromosome for consistent processing
    for chrom in regions:
        regions[chrom].sort()
                    
    total_regions = sum(len(r) for r in regions.values())
    logging.info(f"Read {total_regions} consensus regions with sample information")
    return dict(regions), region_samples

def read_regions_parallel(file_pattern: str, sample_list: List[str], haplotype: str) -> Dict[str, RegionIndex]:
    """Read methylation regions for all samples in parallel."""
    region_indices = {}
    sample_files = []
    
    for sample in sample_list:
        sample_file = file_pattern.format(sample=sample, haplotype=haplotype)
        logging.debug(f"Processing file: {sample_file}")
        if os.path.exists(sample_file):
            sample_files.append((sample_file, sample))
        else:
            logging.warning(f"File not found: {sample_file}")
    
    with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        futures = [executor.submit(process_sample_file, args) for args in sample_files]
        
        for future in tqdm(futures, desc="Reading sample files"):
            try:
                sample_id, regions_by_chrom = future.result()
                if regions_by_chrom:  # Only add if we got valid regions
                    regions_index = RegionIndex()
                    for chrom, regions in regions_by_chrom.items():
                        regions_index.add_regions(chrom, regions)
                    region_indices[sample_id] = regions_index
                else:
                    logging.warning(f"No valid regions found for sample {sample_id}")
            except Exception as e:
                logging.error(f"Error processing sample: {e}")
    
    logging.info("Sorting region indices...")
    for regions_index in tqdm(region_indices.values(), desc="Sorting indices"):
        regions_index.sort_all()
    
    logging.info(f"Successfully loaded regions for {len(region_indices)} samples")
    return region_indices

def process_regions(consensus_regions: Dict[str, List[Tuple[int, int]]], 
                   consensus_samples: Dict[str, Set[str]],
                   unmethylated_regions: Dict[str, RegionIndex],
                   methylated_regions: Dict[str, RegionIndex],
                   min_samples: int) -> List[Tuple]:
    """Process consensus regions against methylated and unmethylated sample data."""
    results = []
    total_regions = sum(len(regions) for regions in consensus_regions.values())
    
    logging.info(f"Processing {total_regions} consensus regions")
    logging.info(f"Loaded {len(methylated_regions)} methylated and {len(unmethylated_regions)} unmethylated sample indices")
    
    skipped_no_unmethylated = 0
    skipped_no_methylated = 0
    skipped_inconsistent = 0
    
    with tqdm(total=total_regions, desc="Processing regions") as pbar:
        for chrom in sorted(consensus_regions.keys(), key=lambda x: CHROM_ORDER.get(x, 999)):
            for start, end in consensus_regions[chrom]:
                region_key = f"{chrom}:{start}-{end}"
                expected_unmethylated = consensus_samples.get(region_key, set())
                
                unmethylated_samples = set()
                methylated_samples = set()
                
                # Check which samples show unmethylation
                for sample_id, region_index in unmethylated_regions.items():
                    if region_index.find_partial_overlap(chrom, start, end):
                        unmethylated_samples.add(sample_id)
                
                # Check which samples have complete methylation coverage
                for sample_id, region_index in methylated_regions.items():
                    if region_index.find_containing_region(chrom, start, end):
                        methylated_samples.add(sample_id)
                
                # Log the methylated samples for the specific region
                if chrom == 'chr16' and start == 3363894 and end == 3363946:
                    logging.debug(f"Region {region_key} methylated samples: {methylated_samples}")
                
                # Verify consistency with consensus
                if not unmethylated_samples.issuperset(expected_unmethylated):
                    skipped_inconsistent += 1
                    logging.debug(f"Region {region_key} missing expected unmethylated samples: "
                                f"{expected_unmethylated - unmethylated_samples}")
                    pbar.update(1)
                    continue
                
                # Track filtering reasons
                if len(unmethylated_samples) < 1:
                    skipped_no_unmethylated += 1
                elif len(methylated_samples) < min_samples:
                    skipped_no_methylated += 1
                
                # Only keep regions where:
                # 1. At least one sample is unmethylated (consensus requirement)
                # 2. Enough samples show complete methylation
                if len(unmethylated_samples) >= 1 and len(methylated_samples) >= min_samples:
                    logging.debug(f"Region {region_key} has {len(methylated_samples)} methylated samples and {len(unmethylated_samples)} unmethylated samples")
                    results.append((
                        chrom,
                        start,
                        end,
                        end - start,
                        len(unmethylated_samples),
                        len(methylated_samples),
                        ','.join(sorted(unmethylated_samples)) if unmethylated_samples else ".",
                        ','.join(sorted(methylated_samples)) if methylated_samples else "."
                    ))
                
                pbar.update(1)
    
    logging.info(f"\nFiltering Summary:")
    logging.info(f"Total input regions: {total_regions}")
    logging.info(f"Regions with no unmethylated samples: {skipped_no_unmethylated}")
    logging.info(f"Regions with insufficient methylated samples (<{min_samples}): {skipped_no_methylated}")
    logging.info(f"Regions with inconsistent unmethylated samples: {skipped_inconsistent}")
    logging.info(f"Final regions after filtering: {len(results)}")
    logging.info(f"Retention rate: {(len(results)/total_regions*100):.2f}%")
    
    return results

def write_results(results: List[Tuple], output_file: str):
    """Write filtered regions to output file."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write("chrom\tstart\tend\tsize\tnum_unmethylated\tnum_methylated\tunmethylated_samples\tmethylated_samples\n")
        for r in results:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\n")

def process_haplotype(args: argparse.Namespace) -> None:
    """Process a single haplotype and generate its output file."""
    # Read sample list
    with open(args.sample_list_file) as f:
        sample_list = [line.strip() for line in f if line.strip()]
    logging.info(f"Read {len(sample_list)} samples")
    
    # Read consensus regions with sample information
    consensus_regions, consensus_samples = read_consensus_regions(args.consensus_regions)
    
    # Create file patterns
    methylated_pattern = os.path.join(args.regions_dir, "{sample}_{haplotype}_M.bed")
    unmethylated_pattern = os.path.join(os.path.dirname(args.regions_dir), f"{args.haplotype}_U", "{sample}_{haplotype}_U.bed")
    
    logging.info(f"Reading methylated regions...")
    methylated_regions = read_regions_parallel(methylated_pattern, sample_list, args.haplotype)
    
    logging.info(f"Reading unmethylated regions...")
    unmethylated_regions = read_regions_parallel(unmethylated_pattern, sample_list, args.haplotype)
    
    # Process regions
    results = process_regions(consensus_regions, consensus_samples, 
                            unmethylated_regions, methylated_regions, args.min_samples)
    
    # Write results
    write_results(results, args.output_file)
    logging.info(f"Completed. Wrote {len(results)} regions to output")

def main():
    parser = argparse.ArgumentParser(description="Filter consensus regions based on methylation coverage.")
    parser.add_argument("--consensus_regions", required=True,
                      help="Path to consensus regions BED file")
    parser.add_argument("--regions_dir", required=True,
                      help="Base directory containing methylation region files")
    parser.add_argument("--output_file", required=True,
                      help="Path to output file")
    parser.add_argument("--sample-list-file", required=True,
                      help="File containing list of sample IDs")
    parser.add_argument("--haplotype", required=True, choices=['H1', 'H2'],
                      help="Haplotype to process (H1 or H2)")
    parser.add_argument("--min-samples", type=int, default=1,
                      help="Minimum number of samples required with complete methylation coverage")
    args = parser.parse_args()
    
    process_haplotype(args)

if __name__ == "__main__":
    main()