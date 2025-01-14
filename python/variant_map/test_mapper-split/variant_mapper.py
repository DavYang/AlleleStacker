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
from typing import Dict, List, Set, Optional, Tuple
from collections import defaultdict
import gzip
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
from pathlib import Path
from tqdm import tqdm
from scipy.stats import fisher_exact  # Import Fisher's exact test

# Import VariantCounts from variant_utils
from variant_utils import VariantCounts

class VariantMethylationMapper:
    def __init__(self, output_prefix: str, haplotype: str, max_workers: Optional[int] = None):
        self.output_prefix = output_prefix
        self.haplotype = haplotype
        self.vcf_files = {}
        self.vcf_samples = {}
        self.sample_map = {}
        self.stats = defaultdict(lambda: defaultdict(int))
        self.max_workers = max_workers or min(mp.cpu_count(), 20)

        # Setup logging
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(f"{output_prefix}.log"),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # Define phased status
        self.phased_types = {
            'small': True,
            'cnv': False,
            'sv': True,
            'tr': False
        }

    def _normalize_sample_name(self, sample: str) -> str:
        """Normalize sample name with special case handling"""
        import re
        return re.sub(r'[_-]', '', sample) 

    def _get_samples(self, sample_str: str) -> Set[str]:
        """Parse and normalize sample string"""
        if pd.isna(sample_str) or sample_str == '.':
            return set()
        return {self._normalize_sample_name(s.strip())
                for s in sample_str.replace(';', ',').split(',') if s.strip()}

    def load_vcfs(self, vcf_files: Dict[str, str]) -> bool:
        """Load and validate VCF files"""
        for vcf_type, path in vcf_files.items():
            if not path or not all(os.path.exists(f) for f in [path, f"{path}.tbi"]):
                self.logger.warning(f"VCF file or index missing for {vcf_type}: {path}")
                continue

            try:
                with gzip.open(path, 'rt') as f:
                    header_found = False
                    for line in f:
                        if line.startswith('#CHROM'):
                            samples = line.strip().split('\t')[9:]
                            self.vcf_samples[vcf_type] = set(samples)
                            self.vcf_files[vcf_type] = path
                            
                            # Build sample mapping
                            for s in samples:
                                normalized_name = self._normalize_sample_name(s)
                                self.sample_map[s] = normalized_name
                                self.sample_map[s.replace('-', '_')] = normalized_name
                            header_found = True
                            break
                    
                    if not header_found:
                        self.logger.error(f"No header found in {vcf_type} VCF")
                        continue
                        
                self.logger.info(f"Loaded {vcf_type} VCF with {len(samples)} samples")
            except Exception as e:
                self.logger.error(f"Error loading {vcf_type} VCF: {e}")
                
        return bool(self.vcf_files)

    def process_variant(self, fields: List[str], vcf_type: str,
                       meth_samples: Set[str], unmeth_samples: Set[str]) -> Optional[Dict]:
        """Process variant with overlap tracking and scoring"""

        genotypes = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
        
        # Initialize counts with all states
        counts = VariantCounts(
            total_meth=len(meth_samples),
            total_unmeth=len(unmeth_samples),
            total_no_data=len(self.vcf_samples[vcf_type] - meth_samples - unmeth_samples)
        )

        meth_vars = []
        unmeth_vars = []
        no_data_vars = []  # Samples with no methylation data, but have variant data
        no_meth_vars = []   # Samples with no methylation data, but have variant data (duplicate for clarity)
        no_var_meth_samples = []  # Methylated samples with no variant data
        no_var_unmeth_samples = [] # Unmethylated samples with no variant data
        no_data_both_samples = [] # Samples with no data for both methylation and variant

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
                        if allele > 0:
                            has_variant = True
                    else:
                        self.logger.warning(f"Ambiguous phasing information for {vcf_type} variant: {fields[2]}")
                except ValueError:
                    self.logger.warning(f"Missing or invalid phasing information for {vcf_type} variant: {fields[2]}")

            if not self.phased_types[vcf_type] or (self.phased_types[vcf_type] and '|' not in gt):
                has_variant = any(int(a) > 0 for a in gt.replace('|', '/').split('/') if a.isdigit())

            if has_variant:
                if sample in meth_samples:
                    meth_vars.append(f"{sample}:{gt}")
                    counts.meth_with_var += 1
                elif sample in unmeth_samples:
                    unmeth_vars.append(f"{sample}:{gt}")
                    counts.unmeth_with_var += 1
                else:
                    no_data_vars.append(f"{sample}:{gt}")
                    no_meth_vars.append(f"{sample}:{gt}")  # Duplicate for clarity
                    counts.no_data_with_var += 1
                    counts.no_meth_with_var += 1

        # Get variant type 
        if vcf_type == 'small':
            if len(fields[3]) == 1 and len(fields[4]) == 1:
                var_type = 'snp'
            elif len(fields[3]) > len(fields[4]):
                var_type = 'deletion'
            else:
                var_type = 'insertion'
        else:
            var_type = vcf_type

        # Calculate association strength metrics
        data_coverage = (counts.total_meth + counts.total_unmeth) / (counts.total_samples - counts.no_data_both)
        meth_association = counts.get_meth_ratio()
        unmeth_association = counts.get_unmeth_ratio()
        
        # Fisher's exact test for significance
        odds_ratio, p_value = fisher_exact([[counts.meth_with_var, counts.total_meth - counts.meth_with_var],
                                          [counts.unmeth_with_var, counts.total_unmeth - counts.unmeth_with_var]])

        # Calculate scores with confidence thresholds and p-value
        if data_coverage >= HIGH_CONFIDENCE:
            confidence_weight = 1.0
        elif data_coverage >= MEDIUM_CONFIDENCE:
            confidence_weight = 0.7
        elif data_coverage >= LOW_CONFIDENCE:
            confidence_weight = 0.4
        else:
            confidence_weight = 0.1

        methylation_score = (confidence_weight * (
                             0.4 * meth_association + 
                             0.6 * (1 - unmeth_association)
                             ) - 0.1 * np.log10(p_value)
                            )

        unmethylation_score = (confidence_weight * (
                               0.4 * unmeth_association + 
                               0.6 * (1 - meth_association)
                               ) - 0.1 * np.log10(p_value)
                              )

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
            'meth_association': meth_association,
            'unmeth_association': unmeth_association,
            'p_value': p_value,
            'meth_samples': meth_vars,
            'unmeth_samples': unmeth_vars,
            'no_data_samples': no_data_vars,  # Samples with no methylation data, but have variant data
            'no_meth_samples': no_meth_vars,   # Samples with no methylation data, but have variant data (duplicate for clarity)
            'no_var_meth_samples': no_var_meth_samples,  # Methylated samples with no variant data
            'no_var_unmeth_samples': no_var_unmeth_samples, # Unmethylated samples with no variant data
            'no_data_both_samples': no_data_both_samples # Samples with no data for both methylation and variant
        }

    def process_region(self, region: pd.Series) -> List[Dict]:
        """Process variants in a region"""
        meth_samples = self._get_samples(region['methylated_samples'])
        unmeth_samples = self._get_samples(region['unmethylated_samples'])

        results = []
        for vcf_type, vcf_path in self.vcf_files.items():
            try:
                cmd = ['tabix', vcf_path, f"{region['chrom']}:{region['start']}-{region['end']}"]
                output = subprocess.check_output(cmd, text=True)

                for line in output.splitlines():
                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        continue

                    variant = self.process_variant(fields, vcf_type, meth_samples, unmeth_samples)
                    if variant:
                        results.append(variant)

            except subprocess.CalledProcessError as e:
                self.logger.warning(f"Error running tabix: {e}")  # More specific error handling
            except Exception as e:
                self.logger.error(f"Error processing {vcf_type} variants: {str(e)}")

        return results

    def run_mapping(self, bed_path: str) -> List[Dict]:
        """Map variants across regions"""
        try:
            df = pd.read_csv(bed_path, sep='\t')
            self.logger.info(f"Processing {len(df)} regions")

            results = []
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
            return results
        
        except FileNotFoundError as e:
            self.logger.error(f"BED file not found: {e}")  # Specific file error
            return []
        except pd.errors.EmptyDataError as e:
            self.logger.error(f"Empty or invalid BED file: {e}")  # Data error
            return []
        except Exception as e:
            self.logger.error(f"Mapping failed: {str(e)}")
            return []

    def write_results(self, variants: List[Dict]) -> bool:
        """Write results with ranking and dual scoring, split by haplotype"""
        try:
            # No sorting here, sorting will be done later

            # Split variants by haplotype
            hap1_variants = [v for v in variants if self.haplotype == 'H1']
            hap2_variants = [v for v in variants if self.haplotype == 'H2']

            # Write variants to separate files for each haplotype
            for hap, hap_variants in [('H1', hap1_variants), ('H2', hap2_variants)]:
                output_file = Path(self.output_prefix).parent / f"{Path(self.output_prefix).stem}_{hap}_scored_variants.tsv"

                with open(output_file, 'w') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow([
                        'chrom', 'start', 'end', 'variant_id', 'type', 'ref', 'alt', 
                        'methylation_score', 'unmethylation_score', 'p_value', 
                        'M+V+', 'M-V+', 'M+V-', 'M-V-', 
                        'X(M-)V+', 'X(V-)M+', 'X(V-)M-', 'X(M-)X(V-)',
                        'meth_samples', 'unmeth_samples', 'no_data_samples',
                        'no_meth_samples', 'no_var_meth_samples', 'no_var_unmeth_samples', 'no_data_both_samples'
                    ])

                    # Sort and write variants for methylation association
                    hap_variants.sort(key=lambda v: (-v['methylation_score'], v['p_value']))  # Sort by methylation score
                    writer.writerow([])  # Add an empty row for separation
                    writer.writerow(["# Methylation-associated Variants"])
                    for var in hap_variants:
                        counts = var['counts']
                        writer.writerow([
                            var['region']['chrom'],
                            var['region']['start'],
                            var['region']['end'],
                            var['variant']['id'],
                            var['variant']['type'],
                            var['variant']['ref'],
                            var['variant']['alt'],
                            var['methylation_score'],  # Use methylation_score here
                            var['unmethylation_score'],  # Include unmethylation score as well
                            var['p_value'],
                            f"{counts.meth_with_var}/{counts.total_meth} ({var['meth_association']:.1%})",
                            f"{counts.unmeth_with_var}/{counts.total_unmeth} ({var['unmeth_association']:.1%})",
                            f"{counts.total_meth - counts.meth_with_var}/{counts.total_meth} ({1 - var['meth_association']:.1%})",
                            f"{counts.total_unmeth - counts.unmeth_with_var}/{counts.total_unmeth} ({1 - var['unmeth_association']:.1%})",
                            f"{counts.no_meth_with_var}/{counts.total_no_data} ({counts.no_meth_with_var / counts.total_no_data:.1%} if counts.total_no_data else '0.0%'}",
                            f"{counts.no_var_meth}/{counts.total_methylated} ({counts.no_var_meth / counts.total_methylated:.1%} if counts.total_methylated else '0.0%'}",
                            f"{counts.no_var_unmeth}/{counts.total_unmethylated} ({counts.no_var_unmeth / counts.total_unmethylated:.1%} if counts.total_unmethylated else '0.0%'}",
                            f"{counts.no_data_both}/{counts.total_samples} ({counts.no_data_both / counts.total_samples:.1%} if counts.total_samples else '0.0%'}",
                            ','.join(var['meth_samples']) or '.',
                            ','.join(var['unmeth_samples']) or '.',
                            ','.join(var['no_data_samples']) or '.',
                            ','.join(var['no_meth_samples']) or '.',
                            ','.join(var['no_var_meth_samples']) or '.',
                            ','.join(var['no_var_unmeth_samples']) or '.',
                            ','.join(var['no_data_both_samples']) or '.'
                        ])

                    # Sort and write variants for unmethylation association
                    hap_variants.sort(key=lambda v: (-v['unmethylation_score'], v['p_value']))  # Sort by unmethylation score
                    writer.writerow([])  # Add an empty row for separation
                    writer.writerow(["# Unmethylation-associated Variants"])
                    for var in hap_variants:
                        counts = var['counts']
                        writer.writerow([
                            var['region']['chrom'],
                            var['region']['start'],
                            var['region']['end'],
                            var['variant']['id'],
                            var['variant']['type'],
                            var['variant']['ref'],
                            var['variant']['alt'],
                            var['methylation_score'],  # Include methylation score as well
                            var['unmethylation_score'],  # Use unmethylation_score here
                            var['p_value'],
                            f"{counts.meth_with_var}/{counts.total_meth} ({var['meth_association']:.1%})",
                            f"{counts.unmeth_with_var}/{counts.total_unmeth} ({var['unmeth_association']:.1%})",
                            f"{counts.total_meth - counts.meth_with_var}/{counts.total_meth} ({1 - var['meth_association']:.1%})",
                            f"{counts.total_unmeth - counts.unmeth_with_var}/{counts.total_unmeth} ({1 - var['unmeth_association']:.1%})",
                            f"{counts.no_meth_with_var}/{counts.total_no_data} ({counts.no_meth_with_var / counts.total_no_data:.1%} if counts.total_no_data else '0.0%'}",
                            f"{counts.no_var_meth}/{counts.total_methylated} ({counts.no_var_meth / counts.total_methylated:.1%} if counts.total_methylated else '0.0%'}",
                            f"{counts.no_var_unmeth}/{counts.total_unmethylated} ({counts.no_var_unmeth / counts.total_unmethylated:.1%} if counts.total_unmethylated else '0.0%'}",
                            f"{counts.no_data_both}/{counts.total_samples} ({counts.no_data_both / counts.total_samples:.1%} if counts.total_samples else '0.0%'}",
                            ','.join(var['meth_samples']) or '.',
                            ','.join(var['unmeth_samples']) or '.',
                            ','.join(var['no_data_samples']) or '.',
                            ','.join(var['no_meth_samples']) or '.',
                            ','.join(var['no_var_meth_samples']) or '.',
                            ','.join(var['no_var_unmeth_samples']) or '.',
                            ','.join(var['no_data_both_samples']) or '.'
                        ])

                self.logger.info(f"Wrote {len(hap_variants)} prioritized variants for {hap} to {output_file}")

            # Write summary statistics
            stats_file = base_path.parent / f"{base_path.stem}_stats.txt"
            with open(stats_file, 'w') as f:
                f.write("Variant Pattern Statistics\n")
                f.write("=========================\n\n")
                
                for vcf_type, patterns in self.stats.items():
                    f.write(f"\n{vcf_type.upper()} Variants:\n")
                    for pattern, count in patterns.items():
                        if pattern != 'total':
                            total = patterns['total']
                            percentage = (count / total * 100) if total > 0 else 0
                            f.write(f"  {pattern}: {count} ({percentage:.1f}%)\n")
            
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
                self.logger.error("No variants found")
                return False

            # Write results
            return self.write_results(variants)

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return False

def main():
    parser = argparse.ArgumentParser(
        description='Map and categorize variants based on methylation patterns',
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