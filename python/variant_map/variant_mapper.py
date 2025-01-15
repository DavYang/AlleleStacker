#!/usr/bin/env python3
"""
variant_methylation_mapper.py - Maps variants to methylation regions with overlap awareness
Handles varied overlap patterns between methylation data and variant calls.
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
from scipy.stats import fisher_exact

# Constants for confidence thresholds
HIGH_CONFIDENCE = 0.8
MEDIUM_CONFIDENCE = 0.6
LOW_CONFIDENCE = 0.4

@dataclass
class VariantCounts:
    """Track sample counts and overlap patterns"""
    total_meth: int = 0
    total_unmeth: int = 0
    total_no_data: int = 0
    meth_with_var: int = 0
    unmeth_with_var: int = 0
    no_meth_with_var: int = 0   # This should match what we access later
    no_data_with_var: int = 0  # Added this missing field
    no_var_meth: int = 0
    no_var_unmeth: int = 0
    no_data_both: int = 0

    @property
    def total_samples(self) -> int:
        return (self.total_meth + self.total_unmeth + self.total_no_data +
                self.no_meth_with_var + self.no_var_meth + self.no_var_unmeth + 
                self.no_data_both)

    @property
    def total_with_variant(self) -> int:
        return self.meth_with_var + self.unmeth_with_var + self.no_data_with_var  # Fixed this

    @property 
    def total_methylated(self) -> int:
        return self.total_meth + self.no_var_meth

    @property
    def total_unmethylated(self) -> int:
        return self.total_unmeth + self.no_var_unmeth

    def get_meth_ratio(self, include_no_variant_data: bool = False) -> float:
        total_meth = self.total_methylated if include_no_variant_data else self.total_meth
        return self.meth_with_var / total_meth if total_meth > 0 else 0.0

    def get_unmeth_ratio(self, include_no_variant_data: bool = False) -> float:
        total_unmeth = self.total_unmethylated if include_no_variant_data else self.total_unmeth
        return self.unmeth_with_var / total_unmeth if total_unmeth > 0 else 0.0


class VariantMethylationMapper:
    def __init__(self, output_prefix: str, haplotype: str, max_workers: Optional[int] = None):
        self.output_prefix = Path(output_prefix)
        self.haplotype = haplotype
        self.vcf_files: Dict[str, str] = {}
        self.vcf_samples: Dict[str, Set[str]] = {}
        self.sample_map: Dict[str, str] = {}
        self.stats = defaultdict(lambda: defaultdict(int))
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
        
        self.phased_types = {
            'small': True,
            'cnv': False,
            'sv': True,
            'tr': False
        }

    def _normalize_sample_name(self, sample: str) -> str:
        """Normalize sample name with special case handling"""
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

    def process_variant(self, fields: List[str], vcf_type: str,
                       meth_samples: Set[str], unmeth_samples: Set[str]) -> Optional[Dict]:
        """Process variant with overlap tracking and scoring"""
        try:
            genotypes = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
            
            counts = VariantCounts(
                total_meth=len(meth_samples),
                total_unmeth=len(unmeth_samples),
                total_no_data=len(self.vcf_samples[vcf_type] - meth_samples - unmeth_samples)
            )

            meth_vars: List[str] = []
            unmeth_vars: List[str] = []
            no_data_vars: List[str] = []
            no_meth_vars: List[str] = []
            no_var_meth_samples: List[str] = []
            no_var_unmeth_samples: List[str] = []
            no_data_both_samples: List[str] = []

            for sample, gt_str in genotypes.items():
                if gt_str in {'./.', '.|.'} or ':' not in gt_str:
                    if sample in meth_samples:
                        no_var_meth_samples.append(sample)
                        counts.no_var_meth += 1
                    elif sample in unmeth_samples:
                        no_var_unmeth_samples.append(sample)
                        counts.no_var_unmeth += 1
                    else:
                        no_data_both_samples.append(sample)
                        counts.no_data_both += 1
                    continue

                gt = gt_str.split(':')[0]
                has_variant = False

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

                if has_variant:
                    if sample in meth_samples:
                        meth_vars.append(f"{sample}:{gt}")
                        counts.meth_with_var += 1
                    elif sample in unmeth_samples:
                        unmeth_vars.append(f"{sample}:{gt}")
                        counts.unmeth_with_var += 1
                    else:
                        no_data_vars.append(f"{sample}:{gt}")
                        no_meth_vars.append(f"{sample}:{gt}")
                        counts.no_data_with_var += 1
                        counts.no_meth_with_var += 1

            # Get variant type
            if vcf_type == 'small':
                var_type = ('snp' if len(fields[3]) == 1 and len(fields[4]) == 1 
                           else 'deletion' if len(fields[3]) > len(fields[4])
                           else 'insertion')
            else:
                var_type = vcf_type

            # Calculate metrics
            data_coverage = ((counts.total_meth + counts.total_unmeth) / 
                           (counts.total_samples - counts.no_data_both)
                           if counts.total_samples > counts.no_data_both else 0.0)
            
            meth_ratio = counts.get_meth_ratio()
            unmeth_ratio = counts.get_unmeth_ratio()

            # Calculate Fisher's exact test
            table = [
                [counts.meth_with_var, counts.total_meth - counts.meth_with_var],
                [counts.unmeth_with_var, counts.total_unmeth - counts.unmeth_with_var]
            ]
            odds_ratio, p_value = fisher_exact(table)

            # Calculate confidence-weighted scores
            confidence_weight = (1.0 if data_coverage >= HIGH_CONFIDENCE else
                               0.7 if data_coverage >= MEDIUM_CONFIDENCE else
                               0.4 if data_coverage >= LOW_CONFIDENCE else 0.1)

            methylation_score = (confidence_weight * 
                               (0.4 * meth_ratio + 0.6 * (1 - unmeth_ratio)) -
                               0.1 * np.log10(p_value if p_value > 0 else 1e-10))

            unmethylation_score = (confidence_weight * 
                                 (0.4 * unmeth_ratio + 0.6 * (1 - meth_ratio)) -
                                 0.1 * np.log10(p_value if p_value > 0 else 1e-10))

            return {
                'region': {
                    'chrom': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[1]) + len(fields[3]) - 1
                },
                'variant': {
                    'id': fields[2] if fields[2] != '.' else f"{var_type}_{fields[0]}_{fields[1]}",
                    'type': var_type,
                    'pos': int(fields[1]),
                    'ref': fields[3],
                    'alt': fields[4]
                },
                'counts': counts,
                'methylation_score': methylation_score,
                'unmethylation_score': unmethylation_score,
                'data_coverage': data_coverage,
                'meth_ratio': meth_ratio,
                'unmeth_ratio': unmeth_ratio,
                'p_value': p_value,
                'meth_samples': meth_vars,
                'unmeth_samples': unmeth_vars,
                'no_data_samples': no_data_vars,
                'no_meth_samples': no_meth_vars,
                'no_var_meth_samples': no_var_meth_samples,
                'no_var_unmeth_samples': no_var_unmeth_samples,
                'no_data_both_samples': no_data_both_samples
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

            self.logger.info(f"Found {len(results)} variants")
            
        except Exception as e:
            self.logger.error(f"Error in mapping: {str(e)}")
            
        return results

    def write_results(self, variants: List[Dict]) -> bool:
        """Write results with methylation scores and patterns"""
        try:
            # Split variants by haplotype
            hap_variants = [v for v in variants if self.haplotype in ('H1', 'H2')]

            # Prepare output file
            output_file = self.output_prefix.parent / f"{self.output_prefix.stem}_scored_variants.tsv"

            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                
                # Write header
                writer.writerow([
                    'chrom', 'start', 'end', 'type', 'ref', 'alt',
                    'methylation_score', 'unmethylation_score', 'p_value',
                    'M+V+', 'M-V+', 'M+V-', 'M-V-',
                    'X(M-)V+', 'X(V-)M+', 'X(V-)M-', 'X(M-)X(V-)',
                    'meth_samples', 'unmeth_samples', 'no_data_samples',
                    'no_meth_samples', 'no_var_meth_samples', 'no_var_unmeth_samples',
                    'no_data_both_samples'
                ])

                # Write methylation-associated variants
                writer.writerow([])
                writer.writerow(["# Methylation-associated Variants"])
                
                # Sort by methylation score
                meth_sorted = sorted(hap_variants, key=lambda v: (-v['methylation_score'], v['p_value']))
                
                for var in meth_sorted:
                    counts = var['counts']
                    meth_ratio = var['meth_ratio']
                    unmeth_ratio = var['unmeth_ratio']
                    
                    writer.writerow([
                        str(var['region']['chrom']),
                        str(var['region']['start']),
                        str(var['region']['end']),
                        str(var['variant']['type']),
                        str(var['variant']['ref']),
                        str(var['variant']['alt']),
                        f"{var['methylation_score']:.3f}",
                        f"{var['unmethylation_score']:.3f}",
                        f"{var['p_value']:.2e}",
                        f"{counts.meth_with_var}/{counts.total_meth} ({meth_ratio:.1%})",
                        f"{counts.unmeth_with_var}/{counts.total_unmeth} ({unmeth_ratio:.1%})",
                        f"{counts.total_meth - counts.meth_with_var}/{counts.total_meth} ({1-meth_ratio:.1%})",
                        f"{counts.total_unmeth - counts.unmeth_with_var}/{counts.total_unmeth} ({1-unmeth_ratio:.1%})",
                        self._format_ratio(counts.no_meth_with_var, counts.total_no_data),
                        self._format_ratio(counts.no_var_meth, counts.total_methylated),
                        self._format_ratio(counts.no_var_unmeth, counts.total_unmethylated),
                        self._format_ratio(counts.no_data_both, counts.total_samples),
                        ','.join(var['meth_samples']) or '.',
                        ','.join(var['unmeth_samples']) or '.',
                        ','.join(var['no_data_samples']) or '.',
                        ','.join(var['no_meth_samples']) or '.',
                        ','.join(var['no_var_meth_samples']) or '.',
                        ','.join(var['no_var_unmeth_samples']) or '.',
                        ','.join(var['no_data_both_samples']) or '.'
                    ])

                # Write unmethylation-associated variants
                writer.writerow([])
                writer.writerow(["# Unmethylation-associated Variants"])
                
                # Sort by unmethylation score
                unmeth_sorted = sorted(hap_variants, key=lambda v: (-v['unmethylation_score'], v['p_value']))
                
                for var in unmeth_sorted:
                    counts = var['counts']
                    meth_ratio = var['meth_ratio']
                    unmeth_ratio = var['unmeth_ratio']
                    
                    writer.writerow([
                        str(var['region']['chrom']),
                        str(var['region']['start']),
                        str(var['region']['end']),
                        str(var['variant']['id']),
                        str(var['variant']['type']),
                        str(var['variant']['ref']),
                        str(var['variant']['alt']),
                        f"{var['methylation_score']:.3f}",
                        f"{var['unmethylation_score']:.3f}",
                        f"{var['p_value']:.2e}",
                        f"{counts.meth_with_var}/{counts.total_meth} ({meth_ratio:.1%})",
                        f"{counts.unmeth_with_var}/{counts.total_unmeth} ({unmeth_ratio:.1%})",
                        f"{counts.total_meth - counts.meth_with_var}/{counts.total_meth} ({1-meth_ratio:.1%})",
                        f"{counts.total_unmeth - counts.unmeth_with_var}/{counts.total_unmeth} ({1-unmeth_ratio:.1%})",
                        self._format_ratio(counts.no_meth_with_var, counts.total_no_data),
                        self._format_ratio(counts.no_var_meth, counts.total_methylated),
                        self._format_ratio(counts.no_var_unmeth, counts.total_unmethylated),
                        self._format_ratio(counts.no_data_both, counts.total_samples),
                        ','.join(var['meth_samples']) or '.',
                        ','.join(var['unmeth_samples']) or '.',
                        ','.join(var['no_data_samples']) or '.',
                        ','.join(var['no_meth_samples']) or '.',
                        ','.join(var['no_var_meth_samples']) or '.',
                        ','.join(var['no_var_unmeth_samples']) or '.',
                        ','.join(var['no_data_both_samples']) or '.'
                    ])

            # Write statistics summary
            stats_file = self.output_prefix.parent / f"{self.output_prefix.stem}_stats.txt"
            with open(stats_file, 'w') as f:
                f.write("Variant Pattern Statistics\n")
                f.write("=========================\n\n")
                
                for vcf_type, patterns in self.stats.items():
                    f.write(f"\n{vcf_type.upper()} Variants:\n")
                    f.write('-' * (len(vcf_type) + 10) + '\n')
                    
                    if 'total' in patterns:
                        total = patterns['total']
                        for pattern, count in patterns.items():
                            if pattern != 'total':
                                percentage = (count / total * 100) if total > 0 else 0
                                f.write(f"{pattern:15} : {count:6d} ({percentage:5.1f}%)\n")

            self.logger.info(f"Wrote {len(hap_variants)} variants to {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"Error writing results: {str(e)}")
            return False

    def _format_ratio(self, numerator: int, denominator: int) -> str:
        """Format ratio with percentage"""
        if denominator == 0:
            return "0/0 (0.0%)"
        return f"{numerator}/{denominator} ({numerator/denominator:.1%})"

    def run(self, bed_path: str) -> bool:
        """Run complete variant mapping pipeline"""
        try:
            # Map variants
            variants = self.run_mapping(bed_path)
            if not variants:
                self.logger.error("No variants found")
                return False

            # Write results
            return self.write_results(variants)

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return False


def main():
    parser = argparse.ArgumentParser(
        description='Map and score variants based on methylation patterns',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example Usage:
  python variant_methylation_mapper.py \\
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