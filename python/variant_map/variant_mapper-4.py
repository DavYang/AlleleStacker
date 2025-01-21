#!/usr/bin/env python3
"""
variant_methylation_mapper.py

Maps variants to methylation regions and scores methylation associations.
Uses an enrichment-based scoring system incorporating confidence and quality factors.
Handles phased/unphased variants from multiple VCF types (SNPs, CNVs, SVs, TRs).
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
    """Track sample counts across all methylation x variant states"""
    # Total samples by methylation state
    total_meth: int = 0
    total_unmeth: int = 0
    total_no_data: int = 0
    
    # Samples with variant
    meth_with_var: List[str] = None  # M+V+
    unmeth_with_var: List[str] = None  # M-V+
    no_meth_with_var: List[str] = None  # X(M)V+
    
    # Samples without variant
    no_var_meth: int = 0        # M+V-
    no_var_unmeth: int = 0      # M-V-
    no_meth_no_var: int = 0     # X(M)V-
    
    # Special cases - incomplete data
    no_var_data_meth: int = 0   # X(V)M+
    no_var_data_unmeth: int = 0 # X(V)M-
    no_data_both: int = 0       # X(M)X(V)
    
    def __post_init__(self):
        # Initialize lists
        self.meth_with_var = []
        self.unmeth_with_var = []
        self.no_meth_with_var = []
    
    @property
    def total_samples(self) -> int:
        """Total number of samples across all states"""
        return (self.total_meth + self.total_unmeth + self.total_no_data +
                len(self.no_meth_with_var) + self.no_var_meth + self.no_var_unmeth +
                self.no_meth_no_var + self.no_data_both)


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

    def calculate_state_enrichment_score(self, counts: VariantCounts, state: str) -> Optional[Dict]:
        """
        Calculate enrichment score for a specific methylation state
        Returns enrichment metrics or None if insufficient data
        
        Args:
            counts: VariantCounts dataclass with sample counts
            state: Either 'methylated' or 'unmethylated'
        """
        # Get relevant counts based on state
        if state == 'methylated':
            state_var = len(counts.meth_with_var)
            state_no_var = counts.no_var_meth
            total_state = state_var + state_no_var
        else:  # unmethylated
            state_var = len(counts.unmeth_with_var)
            state_no_var = counts.no_var_unmeth
            total_state = state_var + state_no_var

        # Check minimum requirements
        MIN_SAMPLES = 10
        if total_state < MIN_SAMPLES:
            return None

        # Calculate state frequency
        state_freq = state_var / total_state if total_state > 0 else 0

        # Calculate background frequency
        total_complete = (len(counts.meth_with_var) + len(counts.unmeth_with_var) + 
                         counts.no_var_meth + counts.no_var_unmeth)
        total_var = len(counts.meth_with_var) + len(counts.unmeth_with_var)
        
        background_freq = total_var / total_complete if total_complete > 0 else 0

        # Avoid division by zero with small epsilon
        epsilon = 1e-10
        
        # Calculate enrichment using log odds ratio
        log_odds = np.log((state_freq + epsilon) / (1 - state_freq + epsilon)) - \
                   np.log((background_freq + epsilon) / (1 - background_freq + epsilon))
        
        # Calculate confidence factor (C)
        confidence = min(1.0, total_state / MIN_SAMPLES)
        
        # Calculate quality factor (Q)
        complete_data = (len(counts.meth_with_var) + len(counts.unmeth_with_var) + 
                        counts.no_var_meth + counts.no_var_unmeth)
        quality = complete_data / counts.total_samples if counts.total_samples > 0 else 0
        
        # Calculate final enrichment score
        enrichment_score = 1 / (1 + np.exp(-log_odds))  # Convert log odds to probability
        final_score = enrichment_score * confidence * quality
        
        return {
            'state': state,
            'enrichment_score': final_score,
            'log_odds': log_odds,
            'state_frequency': state_freq,
            'background_frequency': background_freq,
            'confidence': confidence,
            'quality': quality,
            'complete_samples': complete_data,
            'total_samples': counts.total_samples
        }

    def determine_methylation_association(self, counts: VariantCounts) -> Dict:
        """
        Determine methylation association using updated scoring system
        """
        # Calculate enrichment scores for both states
        meth_score = self.calculate_state_enrichment_score(counts, 'methylated')
        unmeth_score = self.calculate_state_enrichment_score(counts, 'unmethylated')
        
        # Handle cases with insufficient data
        if not meth_score or not unmeth_score:
            return {
                'classification': 'insufficient_data',
                'reason': 'insufficient_samples',
                'methylated_metrics': meth_score,
                'unmethylated_metrics': unmeth_score
            }
            
        # Determine association based on enrichment scores
        STRONG_ASSOC_THRESHOLD = 0.5
        
        if meth_score['enrichment_score'] > unmeth_score['enrichment_score'] * 2:
            classification = 'methylated_associated'
            primary_metrics = meth_score
        elif unmeth_score['enrichment_score'] > meth_score['enrichment_score'] * 2:
            classification = 'unmethylated_associated'
            primary_metrics = unmeth_score
        else:
            classification = 'no_strong_association'
            primary_metrics = None
            
        return {
            'classification': classification,
            'methylated_metrics': meth_score,
            'unmethylated_metrics': unmeth_score,
            'primary_metrics': primary_metrics
        }

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
            
            # Process each sample
            for sample, gt_str in genotypes.items():
                if gt_str in {'./.', '.|.'} or ':' not in gt_str:
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
                        counts.meth_with_var.append(f"{sample}:{gt}")
                    elif sample in unmeth_samples:
                        counts.unmeth_with_var.append(f"{sample}:{gt}")
                    else:
                        counts.no_meth_with_var.append(f"{sample}:{gt}")
                else:
                    if sample in meth_samples:
                        counts.no_var_meth += 1
                    elif sample in unmeth_samples:
                        counts.no_var_unmeth += 1
                    else:
                        counts.no_meth_no_var += 1

            # Get variant type
            if vcf_type == 'small':
                var_type = ('snp' if len(fields[3]) == 1 and len(fields[4]) == 1 
                           else 'deletion' if len(fields[3]) > len(fields[4])
                           else 'insertion')
            else:
                var_type = vcf_type

            # Calculate association using new scoring system
            association = self.determine_methylation_association(counts)
            
            # Only return variants with strong associations unless debugging
            if association['classification'] == 'no_strong_association':
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
                'association': association
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
                            self.excluded_count += len(region_results) == 0
                            pbar.update(1)
                        except Exception as e:
                            self.logger.error(f"Error processing region: {str(e)}")
                            pbar.update(1)

            self.logger.info(f"Found {len(results)} variants with strong associations")
            self.logger.info(f"Excluded {self.excluded_count} variants")
            
        except Exception as e:
            self.logger.error(f"Error in mapping: {str(e)}")
            
        return results

    def write_results(self, variants: List[Dict]) -> bool:
        """Write results with updated association metrics"""
        try:
            # Create output filename
            output_file = f"{self.output_prefix}_variants.tsv"
            
            # Define headers
            headers = [
                'chrom', 'start', 'end', 'variant_id', 'type', 
                'ref', 'alt', 'classification',
                'enrichment_score', 'log_odds', 'confidence', 'quality',
                'state_frequency', 'background_frequency',
                'complete_samples', 'total_samples',
                'meth_var_samples', 'unmeth_var_samples'
            ]
            
            # Write results
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(headers)
                
                for var in variants:
                    assoc = var['association']
                    primary = assoc['primary_metrics']
                    if primary:  # Only write variants with strong associations
                        writer.writerow([
                            var['region']['chrom'],
                            var['region']['start'],
                            var['region']['end'],
                            var['variant']['id'],
                            var['variant']['type'],
                            var['variant']['ref'],
                            var['variant']['alt'],
                            assoc['classification'],
                            f"{primary['enrichment_score']:.3f}",
                            f"{primary['log_odds']:.3f}",
                            f"{primary['confidence']:.3f}",
                            f"{primary['quality']:.3f}",
                            f"{primary['state_frequency']:.3f}",
                            f"{primary['background_frequency']:.3f}",
                            primary['complete_samples'],
                            primary['total_samples'],
                            ','.join(var['counts'].meth_with_var) or '.',
                            ','.join(var['counts'].unmeth_with_var) or '.'
                        ])

            # Write statistics summary
            stats_file = f"{self.output_prefix}_stats.txt"
            with open(stats_file, 'w') as f:
                f.write(f"Variant Methylation Association Analysis\n")
                f.write(f"=====================================\n\n")
                f.write(f"Total variants processed: {len(variants) + self.excluded_count}\n")
                f.write(f"Variants with strong associations: {len(variants)}\n")
                f.write(f"Variants excluded: {self.excluded_count}\n\n")

                # Association distribution
                classifications = defaultdict(int)
                for var in variants:
                    classifications[var['association']['classification']] += 1
                    
                f.write("Association Distribution:\n")
                for cls, count in classifications.items():
                    f.write(f"  {cls}: {count}\n")
                f.write("\n")

                # Score distribution for variants with strong associations
                scores = []
                for var in variants:
                    if var['association']['primary_metrics']:
                        scores.append(var['association']['primary_metrics']['enrichment_score'])
                
                if scores:
                    f.write("Enrichment Score Distribution:\n")
                    f.write(f"  Mean: {np.mean(scores):.3f}\n")
                    f.write(f"  Median: {np.median(scores):.3f}\n")
                    f.write(f"  Std Dev: {np.std(scores):.3f}\n")
                    f.write(f"  Min: {min(scores):.3f}\n")
                    f.write(f"  Max: {max(scores):.3f}\n\n")

                # Quality metrics
                qualities = [var['association']['primary_metrics']['quality'] 
                           for var in variants if var['association']['primary_metrics']]
                if qualities:
                    f.write("Data Quality Distribution:\n")
                    f.write(f"  Mean quality: {np.mean(qualities):.3f}\n")
                    f.write(f"  Median quality: {np.median(qualities):.3f}\n")
                    f.write(f"  Variants with quality >0.8: {sum(1 for q in qualities if q > 0.8)}\n")

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
