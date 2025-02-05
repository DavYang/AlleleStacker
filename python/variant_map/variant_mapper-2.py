#!/usr/bin/env python3
"""
variant_methylation_mapper.py

Maps variants to methylation regions and scores methylation state associations.
Handles all 9 possible methylation x variant states:

States:
- M+V+: Methylated with variant 
- M-V+: Unmethylated with variant
- M+V-: Methylated without variant
- M-V-: Unmethylated without variant 
- X(M)V+: No methylation data, has variant
- X(M)V-: No methylation data, no variant
- X(V)M+: No variant data, methylated
- X(V)M-: No variant data, unmethylated 
- X(M)X(V): No data for either

Scoring prioritizes variants with clear methylation state associations.
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

# Scoring weights
STATE_ENRICHMENT_WEIGHT = 0.5  # Weight for state-specific enrichment
STATE_DEPLETION_WEIGHT = 0.3   # Weight for state-specific depletion
SIGNIFICANCE_WEIGHT = 0.2      # Weight for statistical significance

@dataclass 
class VariantCounts:
    """Track sample counts across all 9 methylation x variant states"""
    # Total samples by methylation state
    total_meth: int = 0
    total_unmeth: int = 0
    total_no_data: int = 0

    # Samples with variant 
    meth_with_var: int = 0      # M+V+
    unmeth_with_var: int = 0    # M-V+
    no_meth_with_var: int = 0   # X(M)V+

    # Samples without variant
    no_var_meth: int = 0        # M+V-
    no_var_unmeth: int = 0      # M-V-
    no_meth_no_var: int = 0     # X(M)V-  [9th state]
    
    # Special cases - incomplete data
    no_var_data_meth: int = 0   # X(V)M+
    no_var_data_unmeth: int = 0 # X(V)M-
    no_data_both: int = 0       # X(M)X(V)

    @property
    def total_samples(self) -> int:
        """Total number of samples across all states"""
        return (self.total_meth + self.total_unmeth + self.total_no_data +
                self.no_meth_with_var + self.no_var_meth + self.no_var_unmeth +
                self.no_meth_no_var + self.no_data_both)

    @property
    def total_with_variant(self) -> int:
        """Samples with confirmed variant presence"""
        return self.meth_with_var + self.unmeth_with_var + self.no_meth_with_var

    @property
    def total_without_variant(self) -> int:
        """Samples with confirmed variant absence"""
        return self.no_var_meth + self.no_var_unmeth + self.no_meth_no_var

    def get_meth_enrichment(self) -> float:
        """Calculate variant enrichment in methylated samples"""
        if self.total_meth == 0:
            return 0.0
        meth_var_ratio = self.meth_with_var / self.total_meth
        unmeth_var_ratio = self.unmeth_with_var / self.total_unmeth if self.total_unmeth > 0 else 0
        return max(0, meth_var_ratio - unmeth_var_ratio)

    def get_unmeth_enrichment(self) -> float:
        """Calculate variant enrichment in unmethylated samples"""
        if self.total_unmeth == 0:
            return 0.0
        unmeth_var_ratio = self.unmeth_with_var / self.total_unmeth
        meth_var_ratio = self.meth_with_var / self.total_meth if self.total_meth > 0 else 0
        return max(0, unmeth_var_ratio - meth_var_ratio)

    def get_data_coverage(self) -> float:
        """Calculate fraction of samples with complete data"""
        samples_with_data = self.total_meth + self.total_unmeth + self.no_meth_with_var + self.no_meth_no_var
        return samples_with_data / self.total_samples if self.total_samples > 0 else 0.0


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
        
        # Define phased vs unphased variant types
        self.phased_types = {
            'small': True,  # SNPs/indels are phased
            'cnv': False,   # CNVs are unphased
            'sv': True,     # SVs are phased
            'tr': False     # TRs are unphased
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

    def process_variant(self, fields: List[str], vcf_type: str,
                       meth_samples: Set[str], unmeth_samples: Set[str]) -> Optional[Dict]:
        """Process variant and classify samples into all 9 states"""
        try:
            genotypes = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
            
            # Initialize counts and sample lists
            counts = VariantCounts(
                total_meth=len(meth_samples),
                total_unmeth=len(unmeth_samples),
                total_no_data=len(self.vcf_samples[vcf_type] - meth_samples - unmeth_samples)
            )

            # Sample lists for each state
            meth_vars: List[str] = []
            unmeth_vars: List[str] = []
            no_data_vars: List[str] = []
            no_meth_vars: List[str] = []
            no_var_meth_samples: List[str] = []
            no_var_unmeth_samples: List[str] = []
            no_meth_no_var_samples: List[str] = []  # X(M)V- samples
            no_data_both_samples: List[str] = []

            for sample, gt_str in genotypes.items():
                if gt_str in {'./.', '.|.'} or ':' not in gt_str:
                    # No variant data
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

                # Check genotype based on variant type
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

                # Classify sample into appropriate state
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
                else:
                    # No variant - track X(M)V- state
                    if sample not in meth_samples and sample not in unmeth_samples:
                        no_meth_no_var_samples.append(f"{sample}:{gt}")
                        counts.no_meth_no_var += 1

            # Calculate confidence based on data coverage
            data_coverage = counts.get_data_coverage()
            confidence_weight = (1.0 if data_coverage >= HIGH_CONFIDENCE else
                               0.7 if data_coverage >= MEDIUM_CONFIDENCE else
                               0.4 if data_coverage >= LOW_CONFIDENCE else 0.1)

            # Calculate enrichment scores
            meth_enrichment = counts.get_meth_enrichment()
            unmeth_enrichment = counts.get_unmeth_enrichment()

            # Calculate Fisher's exact test
            table = [
                [counts.meth_with_var, counts.total_meth - counts.meth_with_var],
                [counts.unmeth_with_var, counts.total_unmeth - counts.unmeth_with_var]
            ]
            _, p_value = fisher_exact(table)
            significance = -np.log10(p_value if p_value > 0 else 1e-10)

            # Calculate methylation association scores
            methylation_score = confidence_weight * (
                STATE_ENRICHMENT_WEIGHT * meth_enrichment +
                STATE_DEPLETION_WEIGHT * (1 - unmeth_enrichment) +
                SIGNIFICANCE_WEIGHT * significance
            )

            unmethylation_score = confidence_weight * (
                STATE_ENRICHMENT_WEIGHT * unmeth_enrichment +
                STATE_DEPLETION_WEIGHT * (1 - meth_enrichment) +
                SIGNIFICANCE_WEIGHT * significance
            )

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
                    'pos': int(fields[1]),
                    'ref': fields[3],
                    'alt': fields[4]
                },
                'counts': counts,
                'methylation_score': methylation_score,
                'unmethylation_score': unmethylation_score,
                'data_coverage': data_coverage,
                'meth_enrichment': meth_enrichment,
                'unmeth_enrichment': unmeth_enrichment,
                'p_value': p_value,
                'meth_samples': meth_vars,
                'unmeth_samples': unmeth_vars,
                'no_data_samples': no_data_vars,
                'no_meth_samples': no_meth_vars,
                'no_meth_no_var_samples': no_meth_no_var_samples,  # NEW: X(M)V- samples
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

    def write_results(self, variants: List[Dict]) -> bool:
        """Write methylation and unmethylation associated variants to separate files"""
        try:
            # Create output filenames
            meth_file = self.output_prefix.parent / f"{self.output_prefix.stem}_methylated_variants.tsv"
            unmeth_file = self.output_prefix.parent / f"{self.output_prefix.stem}_unmethylated_variants.tsv"

            # Define common headers with all 9 states
            headers = [
                'chrom', 'start', 'end', 'type', 'ref', 'alt',
                'methylation_score', 'unmethylation_score', 'p_value',
                'M+V+', 'M-V+', 'M+V-', 'M-V-',
                'X(M)V+', 'X(M)V-',  # Added X(M)V-
                'X(V)M+', 'X(V)M-', 'X(M)X(V)',
                'data_coverage', 'meth_enrichment', 'unmeth_enrichment',
                'meth_samples', 'unmeth_samples', 'no_data_samples',
                'no_meth_samples', 'no_meth_no_var_samples',
                'no_var_meth_samples', 'no_var_unmeth_samples',
                'no_data_both_samples'
            ]

            # Split variants by association
            meth_variants = [v for v in variants if v['methylation_score'] > v['unmethylation_score']]
            unmeth_variants = [v for v in variants if v['unmethylation_score'] >= v['methylation_score']]

            # Sort by scores
            meth_variants.sort(key=lambda x: (-x['methylation_score'], x['p_value']))
            unmeth_variants.sort(key=lambda x: (-x['unmethylation_score'], x['p_value']))

            # Write methylation-associated variants
            with open(meth_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(headers)
                
                for var in meth_variants:
                    writer.writerow(self._format_variant_row(var))

            # Write unmethylation-associated variants
            with open(unmeth_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(headers)
                
                for var in unmeth_variants:
                    writer.writerow(self._format_variant_row(var))

            # Write statistics summary
            self._write_stats_summary(meth_variants, unmeth_variants)
            
            return True

        except Exception as e:
            self.logger.error(f"Error writing results: {str(e)}")
            return False

    def _format_variant_row(self, variant: Dict) -> List[str]:
        """Format variant data for output"""
        counts = variant['counts']
        return [
            str(variant['region']['chrom']),
            str(variant['region']['start']),
            str(variant['region']['end']),
            str(variant['variant']['type']),
            str(variant['variant']['ref']),
            str(variant['variant']['alt']),
            f"{variant['methylation_score']:.3f}",
            f"{variant['unmethylation_score']:.3f}",
            f"{variant['p_value']:.2e}",
            self._format_ratio(counts.meth_with_var, counts.total_meth),
            self._format_ratio(counts.unmeth_with_var, counts.total_unmeth),
            self._format_ratio(counts.no_var_meth, counts.total_meth),
            self._format_ratio(counts.no_var_unmeth, counts.total_unmeth),
            self._format_ratio(counts.no_meth_with_var, counts.total_no_data),
            self._format_ratio(counts.no_meth_no_var, counts.total_no_data),
            self._format_ratio(counts.no_var_data_meth, counts.total_meth),
            self._format_ratio(counts.no_var_data_unmeth, counts.total_unmeth),
            self._format_ratio(counts.no_data_both, counts.total_samples),
            f"{variant['data_coverage']:.3f}",
            f"{variant['meth_enrichment']:.3f}",
            f"{variant['unmeth_enrichment']:.3f}",
            ','.join(variant['meth_samples']) or '.',
            ','.join(variant['unmeth_samples']) or '.',
            ','.join(variant['no_data_samples']) or '.',
            ','.join(variant['no_meth_samples']) or '.',
            ','.join(variant['no_meth_no_var_samples']) or '.',
            ','.join(variant['no_var_meth_samples']) or '.',
            ','.join(variant['no_var_unmeth_samples']) or '.',
            ','.join(variant['no_data_both_samples']) or '.'
        ]

    def _write_stats_summary(self, meth_variants: List[Dict], unmeth_variants: List[Dict]):
        """Write summary statistics"""
        stats_file = self.output_prefix.parent / f"{self.output_prefix.stem}_stats.txt"
        with open(stats_file, 'w') as f:
            f.write(f"Variant Methylation Analysis Summary\n")
            f.write("================================\n\n")
            
            f.write(f"Total variants processed: {len(meth_variants) + len(unmeth_variants)}\n")
            f.write(f"Methylation-associated variants: {len(meth_variants)}\n")
            f.write(f"Unmethylation-associated variants: {len(unmeth_variants)}\n\n")

            # Write variant type breakdown
            f.write("Variant Type Distribution:\n")
            for vcf_type in self.vcf_files:
                meth_count = sum(1 for v in meth_variants if v['variant']['type'] == vcf_type)
                unmeth_count = sum(1 for v in unmeth_variants if v['variant']['type'] == vcf_type)
                f.write(f"{vcf_type}:\n")
                f.write(f"  Methylation-associated: {meth_count}\n")
                f.write(f"  Unmethylation-associated: {unmeth_count}\n")

            # Write score distribution
            f.write("\nScore Distribution:\n")
            meth_scores = [v['methylation_score'] for v in meth_variants]
            unmeth_scores = [v['unmethylation_score'] for v in unmeth_variants]
            
            f.write("Methylation Scores:\n")
            f.write(f"  Mean: {np.mean(meth_scores):.3f}\n")
            f.write(f"  Median: {np.median(meth_scores):.3f}\n")
            f.write(f"  Std: {np.std(meth_scores):.3f}\n")
            
            f.write("Unmethylation Scores:\n")
            f.write(f"  Mean: {np.mean(unmeth_scores):.3f}\n")
            f.write(f"  Median: {np.median(unmeth_scores):.3f}\n")
            f.write(f"  Std: {np.std(unmeth_scores):.3f}\n")

    def _format_ratio(self, numerator: int, denominator: int) -> str:
        """Format ratio with percentage"""
        if denominator == 0:
            return "0/0 (0.0%)"
        percentage = (numerator / denominator) * 100
        return f"{numerator}/{denominator} ({percentage:.1f}%)"

    def run(self, bed_path: str) -> bool:
        """Run complete variant mapping pipeline"""
        try:
            self.logger.info(f"Processing BED file: {bed_path}")
            df = pd.read_csv(bed_path, sep='\t')
            
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = [executor.submit(self.process_region, region) 
                          for _, region in df.iterrows()]
                
                results = []
                with tqdm(total=len(futures), desc="Processing regions") as pbar:
                    for future in futures:
                        try:
                            region_results = future.result()
                            results.extend(region_results)
                            pbar.update(1)
                        except Exception as e:
                            self.logger.error(f"Error processing region: {str(e)}")
                            pbar.update(1)

            if not results:
                self.logger.error("No variants found")
                return False

            return self.write_results(results)

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return False


def main():
    parser = argparse.ArgumentParser(
        description="Map and score variants based on methylation state associations",
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