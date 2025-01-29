#!/usr/bin/env python3
"""
simplified_variant_methylation_mapper.py

Maps variants to methylation regions focusing on clean methylation patterns.
Only keeps variants with perfect methylation status partitioning between
variant carriers and non-carriers.
"""

import argparse
import subprocess
import pandas as pd
import numpy as np
import os
import logging
import csv
from typing import Dict, List, Set, Optional, Tuple
from collections import defaultdict
import gzip
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
from dataclasses import dataclass
from pathlib import Path
from tqdm import tqdm


@dataclass
class VariantCounts:
    """Track sample counts with focus on clean methylation partitioning"""
    # Samples with variant
    meth_with_var: List[str] = None      # M+V+
    unmeth_with_var: List[str] = None    # M-V+
    
    # Samples without variant
    no_var_meth: int = 0        # M+V-
    no_var_unmeth: int = 0      # M-V-
    
    # Missing/incomplete data (tracked but not used in analysis)
    no_geno_meth: int = 0       # M+X
    no_geno_unmeth: int = 0     # M-X
    
    def __post_init__(self):
        # Initialize lists
        self.meth_with_var = []
        self.unmeth_with_var = []
    
    @property
    def total_samples(self) -> int:
        """Total number of samples with complete data"""
        return (len(self.meth_with_var) + len(self.unmeth_with_var) + 
                self.no_var_meth + self.no_var_unmeth)
    
    @property
    def total_variant_carriers(self) -> int:
        """Total number of samples carrying variant"""
        return len(self.meth_with_var) + len(self.unmeth_with_var)
    
    @property
    def total_non_carriers(self) -> int:
        """Total number of samples without variant"""
        return self.no_var_meth + self.no_var_unmeth


class VariantMethylationMapper:
    def __init__(self, output_prefix: str, haplotype: str, max_workers: Optional[int] = None):
        self.output_prefix = Path(output_prefix)
        self.haplotype = haplotype
        self.vcf_files: Dict[str, str] = {}
        self.vcf_samples: Dict[str, Set[str]] = {}
        self.sample_map: Dict[str, str] = {}
        self.max_workers = max_workers or min(mp.cpu_count(), 20)
        
        # Setup logging
        os.makedirs(self.output_prefix.parent, exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(f"{self.output_prefix}.log"),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # Define phased vs unphased variant types
        self.phased_types = {
            'small': True,   # SNPs/indels are phased
            'cnv': False,    # CNVs are unphased
            'sv': True,      # SVs are phased  
            'tr': False      # TRs are unphased
        }

    def _normalize_sample_name(self, sample: str) -> str:
        """Normalize sample name for consistent matching"""
        return sample.replace('-', '').replace('_', '')

    def _get_samples(self, sample_str: str) -> Set[str]:
        """Parse and normalize sample string"""
        if pd.isna(sample_str) or sample_str == '.':
            return set()
        return {self._normalize_sample_name(s.strip())
                for s in sample_str.replace(';', ',').split(',') if s.strip()}

    def load_vcfs(self, vcf_files: Dict[str, str]) -> bool:
        """Load and validate VCF files"""
        for vcf_type, path in vcf_files.items():
            if not path or not os.path.exists(path):
                self.logger.warning(f"VCF file missing for {vcf_type}: {path}")
                continue
                
            if not os.path.exists(f"{path}.tbi"):
                self.logger.warning(f"VCF index missing for {vcf_type}: {path}.tbi")
                continue

            try:
                with gzip.open(path, 'rt') as f:
                    header_found = False
                    for line in f:
                        if line.startswith('#CHROM'):
                            samples = line.strip().split('\t')[9:]
                            self.vcf_samples[vcf_type] = set(samples)
                            self.vcf_files[vcf_type] = path
                            
                            # Build sample name mapping
                            for s in samples:
                                normalized = self._normalize_sample_name(s)
                                self.sample_map[s] = normalized
                                self.sample_map[s.replace('-', '_')] = normalized
                            
                            header_found = True
                            break
                    
                    if not header_found:
                        self.logger.error(f"No header found in {vcf_type} VCF")
                        continue
                        
                self.logger.info(f"Loaded {vcf_type} VCF with {len(samples)} samples")
                
            except Exception as e:
                self.logger.error(f"Error loading {vcf_type} VCF: {str(e)}")
                continue
                
        return bool(self.vcf_files)

    def check_clean_partitioning(self, counts: VariantCounts) -> Optional[str]:
        """
        Check if variant has clean methylation partitioning.
        Returns:
        - 'methylated_associated': All variant carriers are methylated
        - 'unmethylated_associated': All variant carriers are unmethylated
        - None: Mixed or unclear pattern
        """
        # Get counts for variant carriers
        variant_meth = len(counts.meth_with_var)
        variant_unmeth = len(counts.unmeth_with_var)
        
        # Get counts for non-carriers
        no_variant_meth = counts.no_var_meth
        no_variant_unmeth = counts.no_var_unmeth

        # Minimum sample requirements
        MIN_SAMPLES = 5
        if counts.total_variant_carriers < MIN_SAMPLES or counts.total_non_carriers < MIN_SAMPLES:
            return None

        # Check for clean methylation pattern
        if (variant_meth > 0 and variant_unmeth == 0 and 
            no_variant_meth == 0 and no_variant_unmeth > 0):
            return 'methylated_associated'
        
        # Check for clean unmethylation pattern
        if (variant_meth == 0 and variant_unmeth > 0 and
            no_variant_meth > 0 and no_variant_unmeth == 0):
            return 'unmethylated_associated'
        
        return None

    def process_variant(self, fields: List[str], vcf_type: str,
                       meth_samples: Set[str], unmeth_samples: Set[str]) -> Optional[Dict]:
        """Process variant checking for clean methylation patterns"""
        try:
            genotypes = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
            counts = VariantCounts()
            
            # Process each sample's genotype
            for sample, gt_str in genotypes.items():
                if gt_str in {'./.', '.|.'} or ':' not in gt_str:
                    if sample in meth_samples:
                        counts.no_geno_meth += 1
                    elif sample in unmeth_samples:
                        counts.no_geno_unmeth += 1
                    continue

                gt = gt_str.split(':')[0]
                has_variant = False

                # Handle phased vs unphased variants
                if self.phased_types[vcf_type] and '|' in gt:
                    try:
                        alleles = [int(a) for a in gt.split('|')]
                        if len(alleles) == 2:
                            allele = alleles[0 if self.haplotype == 'H1' else 1]
                            has_variant = allele > 0
                    except (ValueError, IndexError):
                        continue
                else:
                    has_variant = any(int(a) > 0 for a in gt.replace('|', '/').split('/') 
                                    if a.isdigit())

                # Assign to appropriate category
                if has_variant:
                    if sample in meth_samples:
                        counts.meth_with_var.append(f"{sample}:{gt}")
                    elif sample in unmeth_samples:
                        counts.unmeth_with_var.append(f"{sample}:{gt}")
                else:
                    if sample in meth_samples:
                        counts.no_var_meth += 1
                    elif sample in unmeth_samples:
                        counts.no_var_unmeth += 1

            # Check methylation pattern
            classification = self.check_clean_partitioning(counts)
            if not classification:
                return None  # Skip variants without clean patterns

            # Get variant type
            if vcf_type == 'small':
                var_type = ('snp' if len(fields[3]) == 1 and len(fields[4]) == 1 
                           else 'deletion' if len(fields[3]) > len(fields[4])
                           else 'insertion')
            else:
                var_type = vcf_type
                
            return {
                'region': {
                    'chrom': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[1]) + len(fields[3]) - 1
                },
                'variant': {
                    'id': fields[2] if fields[2] != '.' else f"{var_type}_{fields[0]}_{fields[1]}",
                    'type': var_type,
                    'ref': fields[3],
                    'alt': fields[4]
                },
                'classification': classification,
                'metrics': {
                    'variant_carriers': counts.total_variant_carriers,
                    'non_carriers': counts.total_non_carriers,
                    'total_samples': counts.total_samples,
                    'missing_data': counts.no_geno_meth + counts.no_geno_unmeth
                },
                'samples': {
                    'meth_with_var': counts.meth_with_var,
                    'unmeth_with_var': counts.unmeth_with_var,
                }
            }
            
        except Exception as e:
            self.logger.error(f"Error processing variant: {str(e)}")
            return None

    def process_region(self, region: pd.Series) -> List[Dict]:
        """Process variants in a region"""
        results: List[Dict] = []
        
        try:
            meth_samples = self._get_samples(region['methylated_samples'])
            unmeth_samples = self._get_samples(region['unmethylated_samples'])

            for vcf_type, vcf_path in self.vcf_files.items():
                cmd = ['tabix', vcf_path, f"{region['chrom']}:{region['start']}-{region['end']}"]
                
                try:
                    output = subprocess.check_output(cmd, text=True)
                except subprocess.CalledProcessError:
                    continue
                except Exception as e:
                    self.logger.error(f"Error running tabix: {str(e)}")
                    continue

                for line in output.splitlines():
                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        continue

                    variant = self.process_variant(fields, vcf_type, meth_samples, unmeth_samples)
                    if variant:
                        results.append(variant)

        except Exception as e:
            self.logger.error(f"Error processing region: {str(e)}")
            
        return results

    def run_mapping(self, bed_path: str) -> List[Dict]:
        """Map variants across regions"""
        results: List[Dict] = []
        
        try:
            df = pd.read_csv(bed_path, sep='\t')
            self.logger.info(f"Processing {len(df)} regions")

            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = [executor.submit(self.process_region, region) 
                          for _, region in df.iterrows()]
                
                with tqdm(total=len(futures), desc="Processing regions") as pbar:
                    for future in futures:
                        try:
                            region_results = future.result()
                            results.extend(region_results)
                            pbar.update(1)
                        except Exception as e:
                            self.logger.error(f"Error processing region: {str(e)}")
                            pbar.update(1)

            self.logger.info(f"Found {len(results)} variants with clean methylation patterns")
            
        except Exception as e:
            self.logger.error(f"Error in mapping: {str(e)}")
            
        return results

    def write_results(self, variants: List[Dict]) -> bool:
        """Write results focusing on clean methylation patterns"""
        try:
            # Create output filename
            output_file = f"{self.output_prefix}_variants.tsv"
            
            # Define headers
            headers = [
                'chrom', 'start', 'end', 'variant_id', 'type', 
                'ref', 'alt', 'methylation_pattern',
                'variant_carriers', 'non_carriers', 'total_samples', 'missing_data',
                'meth_var_samples', 'unmeth_var_samples'
            ]
            
            # Write results
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(headers)
                
                for var in variants:
                    writer.writerow([
                        var['region']['chrom'],
                        var['region']['start'],
                        var['region']['end'],
                        var['variant']['id'],
                        var['variant']['type'],
                        var['variant']['ref'],
                        var['variant']['alt'],
                        var['classification'],
                        var['metrics']['variant_carriers'],
                        var['metrics']['non_carriers'],
                        var['metrics']['total_samples'],
                        var['metrics']['missing_data'],
                        ','.join(var['samples']['meth_with_var']) or '.',
                        ','.join(var['samples']['unmeth_with_var']) or '.'
                    ])

            # Write statistics summary
            stats_file = f"{self.output_prefix}_stats.txt"
            with open(stats_file, 'w') as f:
                f.write(f"Clean Methylation Pattern Analysis\n")
                # Continue writing statistics file
                f.write(f"================================\n\n")
                f.write(f"Analysis Summary:\n")
                f.write(f"----------------\n")
                
                # Pattern distribution
                pattern_counts = defaultdict(int)
                variant_types = defaultdict(int)
                sample_counts = []
                
                for var in variants:
                    pattern_counts[var['classification']] += 1
                    variant_types[var['variant']['type']] += 1
                    sample_counts.append(var['metrics']['total_samples'])
                
                f.write("\nMethylation Patterns:\n")
                for pattern, count in pattern_counts.items():
                    f.write(f"  {pattern}: {count}\n")
                
                f.write("\nVariant Types:\n")
                for vtype, count in variant_types.items():
                    f.write(f"  {vtype}: {count}\n")
                
                if sample_counts:
                    f.write("\nSample Coverage:\n")
                    f.write(f"  Mean samples per variant: {np.mean(sample_counts):.1f}\n")
                    f.write(f"  Median samples per variant: {np.median(sample_counts):.1f}\n")
                    f.write(f"  Min samples: {min(sample_counts)}\n")
                    f.write(f"  Max samples: {max(sample_counts)}\n")

            self.logger.info(f"Results written to: {output_file}")
            self.logger.info(f"Statistics written to: {stats_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing results: {str(e)}")
            return False

    def run(self, bed_path: str) -> bool:
        """Run complete variant mapping pipeline"""
        try:
            # Map variants
            variants = self.run_mapping(bed_path)
            if not variants:
                self.logger.error("No variants with clean methylation patterns found")
                return False

            # Write results
            return self.write_results(variants)

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return False


def main():
    parser = argparse.ArgumentParser(
        description="""
Map variants to methylation regions focusing on clean methylation patterns.
Only reports variants where methylation status perfectly partitions between
variant carriers and non-carriers.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example Usage:
  python simplified_variant_methylation_mapper.py \\
    --bed regions.bed \\
    --haplotype H1 \\
    --output-prefix results/variants \\
    --small-vcf variants.vcf.gz \\
    --cnv-vcf cnvs.vcf.gz \\
    --sv-vcf svs.vcf.gz \\
    --tr-vcf repeats.vcf.gz
        """
    )
    
    # Required arguments
    parser.add_argument('--bed', required=True,
                       help='Input BED file with methylation regions')
    parser.add_argument('--haplotype', required=True, choices=['H1', 'H2'],
                       help='Haplotype to process')
    parser.add_argument('--output-prefix', required=True,
                       help='Prefix for output files')
    
    # Optional variant files
    parser.add_argument('--small-vcf',
                       help='Small variants (SNPs/indels) VCF file')
    parser.add_argument('--cnv-vcf',
                       help='Copy Number Variants (CNV) VCF file')
    parser.add_argument('--sv-vcf',
                       help='Structural Variants (SV) VCF file')
    parser.add_argument('--tr-vcf',
                       help='Tandem Repeats (TR) VCF file')
    
    # Optional processing parameters
    parser.add_argument('--threads', type=int,
                       help='Number of processing threads (default: auto)')
    
    args = parser.parse_args()
    
    # Initialize mapper
    mapper = VariantMethylationMapper(
        output_prefix=args.output_prefix,
        haplotype=args.haplotype,
        max_workers=args.threads
    )
    
    # Prepare VCF files dictionary
    vcf_files = {
        'small': args.small_vcf,
        'cnv': args.cnv_vcf,
        'sv': args.sv_vcf,
        'tr': args.tr_vcf
    }
    
    # Load VCF files
    if not mapper.load_vcfs(vcf_files):
        mapper.logger.error("No valid VCF files provided")
        exit(1)
    
    # Run pipeline
    if not mapper.run(args.bed):
        mapper.logger.error("Variant mapping failed")
        exit(1)
    
    mapper.logger.info("Variant mapping completed successfully")


if __name__ == '__main__':
    main()
