#!/usr/bin/env python3
"""
variant_methylation_mapper.py

Maps variants to methylation regions and scores methylation associations.
Handles phased/unphased variants from multiple VCF types (SNPs, CNVs, SVs, TRs) 
with complete state tracking across 9 methylation x variant states.
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
    no_meth_no_var: int = 0     # X(M)V-
    
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
    
    def get_data_coverage(self) -> float:
        """Calculate fraction of samples with complete data"""
        samples_with_data = self.total_meth + self.total_unmeth + self.no_meth_with_var + self.no_meth_no_var
        return samples_with_data / self.total_samples if self.total_samples > 0 else 0.0
    
    def get_meth_ratio(self) -> float:
        """Calculate variant presence ratio in methylated samples"""
        return self.meth_with_var / self.total_meth if self.total_meth > 0 else 0.0
    
    def get_unmeth_ratio(self) -> float:
        """Calculate variant presence ratio in unmethylated samples"""
        return self.unmeth_with_var / self.total_unmeth if self.total_unmeth > 0 else 0.0


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
        
    def calculate_methylation_score(self, counts: VariantCounts) -> float:
        """
        Calculate normalized methylation score (0-1)
        0.0 = Strongly unmethylated
        0.5 = No clear association
        1.0 = Strongly methylated
        """
        if counts.meth_with_var + counts.unmeth_with_var == 0:
            return None

        # Calculate variant meth rate
        variant_meth_rate = counts.meth_with_var / (counts.meth_with_var + counts.unmeth_with_var)
        
        # Calculate background meth rate 
        if counts.no_var_meth + counts.no_var_unmeth > 0:
            background_meth_rate = counts.no_var_meth / (counts.no_var_meth + counts.no_var_unmeth)
        else:
            background_meth_rate = 0.0

        # Calculate consistency
        total_with_variant = counts.meth_with_var + counts.unmeth_with_var
        consistency = (1 - np.std([1] * counts.meth_with_var + [0] * counts.unmeth_with_var) / 
                      np.sqrt(total_with_variant)) if total_with_variant > 0 else 0.0

        # Calculate sample confidence
        sample_confidence = min(1.0, (counts.meth_with_var + counts.unmeth_with_var) / 10.0)

        # Calculate data coverage confidence
        data_coverage = counts.get_data_coverage()
        coverage_confidence = min(1.0, data_coverage * 2)  # Full confidence at 50% coverage

        # Calculate enrichment component (-1 to 1)
        enrichment = variant_meth_rate - background_meth_rate

        # Combine components into final score (removed significance component)
        raw_score = (
            0.5 * enrichment +          # Increased weight for enrichment
            0.3 * consistency +         # Pattern consistency
            0.2 * coverage_confidence   # Data coverage 
        )

        # Scale -1,1 range to 0,1 range
        score = (raw_score + 1) / 2

        # Apply confidence 
        final_score = score * sample_confidence

        return max(0.0, min(1.0, final_score))

    def process_variant(self, fields: List[str], vcf_type: str,
                       meth_samples: Set[str], unmeth_samples: Set[str]) -> Optional[Dict]:
        """Process variant with complete state tracking"""
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
            no_meth_no_var_samples: List[str] = []
            no_data_both_samples: List[str] = []
            
            # Process each sample
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
                        counts.no_meth_with_var += 1
                else:
                    if sample not in meth_samples and sample not in unmeth_samples:
                        no_meth_no_var_samples.append(f"{sample}:{gt}")
                        counts.no_meth_no_var += 1

            # Get variant type
            if vcf_type == 'small':
                var_type = ('snp' if len(fields[3]) == 1 and len(fields[4]) == 1 
                           else 'deletion' if len(fields[3]) > len(fields[4])
                           else 'insertion')
            else:
                var_type = vcf_type

            # Calculate methylation score (no p-value used)
            methylation_score = self.calculate_methylation_score(counts)
            if methylation_score is None:
                return None

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
                'samples': {
                    'meth_with_var': meth_vars,
                    'unmeth_with_var': unmeth_vars,
                    'no_data_with_var': no_data_vars,
                    'no_meth_with_var': no_meth_vars,
                    'no_var_meth': no_var_meth_samples,
                    'no_var_unmeth': no_var_unmeth_samples,
                    'no_meth_no_var': no_meth_no_var_samples,
                    'no_data_both': no_data_both_samples
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
        self.excluded_count = 0  # Track excluded variants
        
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
                            self.excluded_count += len(region_results) == 0  # Increment if no variants passed
                            pbar.update(1)
                        except Exception as e:
                            self.logger.error(f"Error processing region: {str(e)}")
                            pbar.update(1)

            self.logger.info(f"Found {len(results)} variants")
            self.logger.info(f"Excluded {self.excluded_count} variants")
            
        except Exception as e:
            self.logger.error(f"Error in mapping: {str(e)}")
            
        return results

    def write_results(self, variants: List[Dict]) -> bool:
        """Write results with methylation scores and patterns"""
        try:
            # Create output filename
            output_file = f"{self.output_prefix}_variants.tsv"

            # Define headers (removed p-value)
            headers = [
                'chrom', 'start', 'end', 'variant_id', 'type', 
                'ref', 'alt', 'methylation_score', 
                'M+V+', 'M-V+', 'M+V-', 'M-V-', 'X(M)V+', 'X(M)V-', 
                'X(V)M+', 'X(V)M-', 'X(M)X(V)',
                'meth_samples', 'unmeth_samples', 'no_data_samples'
            ]

            # Sort variants by methylation score
            variants.sort(key=lambda x: x['methylation_score'], reverse=True)

            # Write results
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(headers)
                
                for var in variants:
                    counts = var['counts']
                    writer.writerow([
                        var['region']['chrom'],
                        var['region']['start'],
                        var['region']['end'],
                        var['variant']['id'],
                        var['variant']['type'],
                        var['variant']['ref'],
                        var['variant']['alt'],
                        f"{var['methylation_score']:.3f}",  # Removed p-value output
                        f"{counts.meth_with_var}/{counts.total_meth}",
                        f"{counts.unmeth_with_var}/{counts.total_unmeth}",
                        f"{counts.no_var_meth}/{counts.total_meth}",
                        f"{counts.no_var_unmeth}/{counts.total_unmeth}",
                        f"{counts.no_meth_with_var}/{counts.total_no_data}",
                        f"{counts.no_meth_no_var}/{counts.total_no_data}",
                        f"{counts.no_var_data_meth}/{counts.total_meth}",
                        f"{counts.no_var_data_unmeth}/{counts.total_unmeth}",
                        f"{counts.no_data_both}/{counts.total_no_data}",
                        ','.join(var['samples']['meth_with_var']) or '.',
                        ','.join(var['samples']['unmeth_with_var']) or '.',
                        ','.join(var['samples']['no_data_with_var']) or '.'
                    ])

            # Write statistics summary
            stats_file = f"{self.output_prefix}_stats.txt"
            with open(stats_file, 'w') as f:
                f.write(f"Variant Methylation Analysis Summary\n")
                f.write(f"=================================\n\n")
                f.write(f"Total variants analyzed: {len(variants) + self.excluded_count}\n")
                f.write(f"Variants passing filters: {len(variants)}\n")
                f.write(f"Variants excluded: {self.excluded_count}\n\n")

                # Overall score distribution
                scores = [v['methylation_score'] for v in variants]
                f.write("Score Distribution:\n")
                f.write(f"  Strongly methylated (>0.7): {sum(1 for s in scores if s > 0.7)}\n")
                f.write(f"  Moderately methylated (0.6-0.7): {sum(1 for s in scores if 0.6 <= s <= 0.7)}\n")
                f.write(f"  No clear association (0.4-0.6): {sum(1 for s in scores if 0.4 <= s <= 0.6)}\n")
                f.write(f"  Moderately unmethylated (0.3-0.4): {sum(1 for s in scores if 0.3 <= s <= 0.4)}\n")
                f.write(f"  Strongly unmethylated (<0.3): {sum(1 for s in scores if s < 0.3)}\n\n")

                # By variant type
                f.write("Variants by Type:\n")
                type_counts = defaultdict(int)
                for var in variants:
                    type_counts[var['variant']['type']] += 1
                for vtype, count in sorted(type_counts.items()):
                    f.write(f"  {vtype}: {count}\n")

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
                self.logger.error("No variants found")
                return False

            # Write results
            return self.write_results(variants)

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return False


def main():
    parser = argparse.ArgumentParser(
        description='Map variants to methylation regions and score methylation associations',
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