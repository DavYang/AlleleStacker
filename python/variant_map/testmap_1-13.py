#!/usr/bin/env python3
"""
Variant Methylation Mapper and Scorer

Comprehensive pipeline to:
1. Map variants across methylation regions
2. Score and prioritize variants based on methylation associations
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
import csv
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
from tqdm import tqdm
from dataclasses import dataclass
from pathlib import Path

@dataclass
class VariantCounts:
    """Container for variant distribution counts"""
    methylated_total: int
    unmethylated_total: int
    no_data_total: int
    methylated_with_var: int
    unmethylated_with_var: int
    no_data_with_var: int
    
    @property
    def total_samples(self) -> int:
        return self.methylated_total + self.unmethylated_total + self.no_data_total
    
    @property
    def total_with_variant(self) -> int:
        return self.methylated_with_var + self.unmethylated_with_var + self.no_data_with_var
    
    @property
    def methylated_ratio(self) -> float:
        return self.methylated_with_var / self.methylated_total if self.methylated_total > 0 else 0
    
    @property
    def unmethylated_ratio(self) -> float:
        return self.unmethylated_with_var / self.unmethylated_total if self.unmethylated_total > 0 else 0

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
        
        # Phased status for variant types
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
                            break
                self.logger.info(f"Loaded {vcf_type} VCF: {len(samples)} samples")
            except Exception as e:
                self.logger.error(f"Error loading {vcf_type} VCF: {e}")
                
        return bool(self.vcf_files)

    def process_variant(self, fields: List[str], vcf_type: str,
                       meth_samples: Set[str], unmeth_samples: Set[str]) -> Optional[Dict]:
        """Process variant checking haplotype phase and categorizing samples"""
        genotypes = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
        meth_vars, unmeth_vars, no_data_vars = [], [], []
        meth_samples_with_var = set()  # Track methylated samples with variant
        unmeth_samples_with_var = set()  # Track unmethylated samples with variant

        for sample, gt_str in genotypes.items():
            if gt_str in {'./.', '.|.'} or ':' not in gt_str:
                no_data_vars.append(sample)
                continue
                
            gt = gt_str.split(':')[0]

            # Check for missing or ambiguous phasing information
            if self.phased_types[vcf_type] and '|' in gt:
                try:
                    # Attempt to split and convert to integers
                    alleles = [int(a) for a in gt.split('|')]
                    if len(alleles) == 2:
                        allele = alleles[0 if self.haplotype == 'H1' else 1]
                        if allele > 0:
                            if sample in meth_samples:
                                meth_vars.append(sample)
                                meth_samples_with_var.add(sample)
                            elif sample in unmeth_samples:
                                unmeth_vars.append(sample)
                                unmeth_samples_with_var.add(sample)
                    else:
                        # Handle ambiguous phasing (e.g., more than 2 alleles) as unphased
                        self.logger.warning(f"Ambiguous phasing information for {vcf_type} variant: {fields[2]}")
                        # Fallback to unphased handling below
                except ValueError:
                    # Handle missing or non-integer alleles as unphased
                    self.logger.warning(f"Missing or invalid phasing information for {vcf_type} variant: {fields[2]}")
                    # Fallback to unphased handling below

            if not self.phased_types[vcf_type] or (self.phased_types[vcf_type] and '|' not in gt):  
                # Unphased - check for any alt allele
                if any(int(a) > 0 for a in gt.replace('|', '/').split('/') if a.isdigit()):
                    if sample in meth_samples:
                        meth_vars.append(sample)
                        meth_samples_with_var.add(sample)
                    elif sample in unmeth_samples:
                        unmeth_vars.append(sample)
                        unmeth_samples_with_var.add(sample)

        if not (meth_vars or unmeth_vars):
            return None

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

        # Update variant info
        variant = {
            'id': fields[2] if fields[2] != '.' else f"{var_type}_{fields[0]}_{fields[1]}",
            'type': var_type,
            'pos': int(fields[1]),
            'ref': fields[3],
            'alt': fields[4],
            'meth_samples': meth_samples_with_var,  # Add the set of samples to the variant info
            'unmeth_samples': unmeth_samples_with_var  # Add the set of samples to the variant info
        }

        # Prepare sample count data
        result = {
            'region': {
                'chrom': fields[0],
                'start': int(fields[1]),
                'end': int(fields[1]) + len(fields[3]) - 1
            },
            'variant': variant,
            'total_meth': len(meth_samples),
            'total_unmeth': len(unmeth_samples),
            'total_no_data': len(no_data_vars),
            'meth_with_var': len(meth_vars),
            'unmeth_with_var': len(unmeth_vars),
            'no_data_with_var': len(no_data_vars)
        }

        return result


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

            except subprocess.CalledProcessError:
                pass
            except Exception as e:
                self.logger.error(f"Error processing {vcf_type} variants: {str(e)}")

        return results

    def run_mapping(self, bed_path: str) -> List[Dict]:
        """Run variant mapping on BED file"""
        try:
            df = pd.read_csv(bed_path, sep='\t')
            self.logger.info(f"Processing {len(df)} regions")

            results = []
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = [executor.submit(self.process_region, region) 
                          for _, region in df.iterrows()]
                
                for future in tqdm(futures, desc="Processing regions"):
                    try:
                        region_results = future.result()
                        results.extend(region_results)
                    except Exception as e:
                        self.logger.error(f"Error processing region: {str(e)}")

            self.logger.info(f"Found {len(results)} potential variants")
            return results
        
        except Exception as e:
            self.logger.error(f"Mapping failed: {str(e)}")
            return []

    def get_effect_category(self, counts: VariantCounts) -> str:
        """Categorize effect size based on ratio differences"""
        ratio_diff = abs(counts.methylated_ratio - counts.unmethylated_ratio)
        
        # Using constants for thresholds
        STRONG_EFFECT = 0.6
        MODERATE_EFFECT = 0.2
        
        if ratio_diff > STRONG_EFFECT:
            return "Strong"
        elif ratio_diff > MODERATE_EFFECT:
            return "Moderate"
        else:
            return "Weak"

    def get_confidence_rating(self, counts: VariantCounts) -> str:
        """Rate confidence based on data completeness"""
        total_with_var = counts.total_with_variant
        if total_with_var < 3:
            return "Low (few samples)"
            
        with_data = counts.methylated_with_var + counts.unmethylated_with_var
        data_ratio = with_data / total_with_var
        
        # Using constants for thresholds
        HIGH_CONFIDENCE = 0.8
        MEDIUM_CONFIDENCE = 0.6
        
        if data_ratio > HIGH_CONFIDENCE:
            return "High"
        elif data_ratio > MEDIUM_CONFIDENCE:
            return "Medium"
        else:
            return "Low"

    def get_tier(self, effect: str, confidence: str) -> int:
        """Determine priority tier based on effect and confidence"""
        if effect == "Strong" and confidence == "High":
            return 1
        elif (effect == "Strong" and confidence == "Medium") or \
             (effect == "Moderate" and confidence == "High"):
            return 2
        elif effect == "Moderate" and confidence == "Medium":
            return 3
        else:
            return 4

    def score_variants(self, mapped_variants: List[Dict]) -> Tuple[List[Dict], List[Dict], List[Dict]]:
        """Score mapped variants based on methylation associations"""
        meth_associated, unmeth_associated, rare_variants = [], [], []  # Add rare_variants list

        for variant in mapped_variants:
            counts = VariantCounts(
                methylated_total=variant['total_meth'],
                unmethylated_total=variant['total_unmeth'],
                no_data_total=variant['total_no_data'],
                methylated_with_var=variant['meth_with_var'],
                unmethylated_with_var=variant['unmeth_with_var'],
                no_data_with_var=variant['no_data_with_var']
            )
            
            # Calculate effect and confidence
            effect = self.get_effect_category(counts)
            confidence = self.get_confidence_rating(counts)
            tier = self.get_tier(effect, confidence)
            
            # Determine bias direction
            bias_direction = 'methylated' if counts.methylated_ratio > counts.unmethylated_ratio else 'unmethylated'
            
            # Prepare scored variant
            scored_variant = {
                # Region info
                'region_chr': variant['region']['chrom'],
                'region_start': variant['region']['start'],
                'region_end': variant['region']['end'],
                'total_samples': counts.total_samples,
                
                # Variant info
                'variant_id': variant['variant']['id'],
                'variant_type': variant['variant']['type'],
                'variant_pos': variant['variant']['pos'],
                'ref': variant['variant']['ref'],
                'alt': variant['variant']['alt'],
                
                # Score info
                'effect': effect,
                'confidence': confidence,
                'tier': tier,
                'bias_direction': bias_direction,
                
                # Ratios
                'meth_ratio': f"{variant['meth_with_var']}/{variant['total_meth']} ({counts.methylated_ratio:.1%})",
                'unmeth_ratio': f"{variant['unmeth_with_var']}/{variant['total_unmeth']} ({counts.unmethylated_ratio:.1%})",
                'no_data_ratio': f"{variant['no_data_with_var']}/{variant['total_no_data']} ({variant['no_data_with_var']/variant['total_no_data']:.1%})",
                
                # Samples
                'meth_samples': ','.join(variant['meth_samples']),  # Add the sample list to the output
                'unmeth_samples': ','.join(variant['unmeth_samples']),  # Add the sample list to the output
                
                # For sorting
                'effect_strength': abs(counts.methylated_ratio - counts.unmethylated_ratio)
            }

            # Check for rare variants with limited data
            if counts.total_with_variant < 3 or (counts.methylated_with_var + counts.unmethylated_with_var) < 2:
                rare_variants.append(scored_variant)  # Add to rare_variants list
            else:
                # Categorize by bias direction as before
                if bias_direction == 'methylated':
                    meth_associated.append(scored_variant)
                else:
                    unmeth_associated.append(scored_variant)
        
        # Sort variants by tier and effect strength
        for variant_list in [meth_associated, unmeth_associated]:
            variant_list.sort(key=lambda x: (x['tier'], -x['effect_strength']))
        
        return meth_associated, unmeth_associated, rare_variants  # Return the rare_variants list

    def write_prioritized_variants(self, meth_variants: List[Dict], unmeth_variants: List[Dict], rare_variants: List[Dict]):
        """Write prioritized variants to separate files with tier organization"""
        columns = [
            'region_chr', 'region_start', 'region_end',
            'variant_id', 'variant_type', 'variant_pos', 
            'ref', 'alt', 'effect', 'confidence', 'tier',
            'meth_ratio', 'unmeth_ratio', 'no_data_ratio',
            'meth_samples', 'unmeth_samples'  # Add the sample columns to the output
        ]
        
        def write_tiered_variants(variants: List[Dict], filename: str):
            with open(filename, 'w') as f:
                # Write header
                f.write("# Variant Prioritization Results\n")
                f.write("# =============================\n\n")
                
                # Refactor to reduce redundancy
                for tier in range(1, 5):
                    tier_vars = [v for v in variants if v['tier'] == tier]
                    if tier_vars:
                        f.write(f"\n## Tier {tier} Variants (n={len(tier_vars)})\n")
                        f.write(pd.DataFrame(tier_vars)[columns].to_string(index=False))
                        f.write("\n")

        # Write rare variants to a separate file
        rare_output = f"{self.output_prefix}_rare_variants.tsv"
        with open(rare_output, 'w') as f:
            f.write("# Rare Variants with Limited Data\n")
            f.write("# ==============================\n\n")
            if rare_variants:
                f.write(pd.DataFrame(rare_variants)[columns].to_string(index=False))
            f.write("\n")

        self.logger.info(f"Wrote {len(rare_variants)} rare variants to {rare_output}")
        
        # Write both files
        meth_output = f"{self.output_prefix}_methylated_variants.tsv"
        unmeth_output = f"{self.output_prefix}_unmethylated_variants.tsv"
        
        write_tiered_variants(meth_variants, meth_output)
        write_tiered_variants(unmeth_variants, unmeth_output)
        
        self.logger.info(f"Wrote {len(meth_variants)} methylation-associated variants to {meth_output}")
        self.logger.info(f"Wrote {len(unmeth_variants)} unmethylation-associated variants to {unmeth_output}")

    def run(self, bed_path: str) -> bool:
        """
        Complete workflow: map variants, score them, and write prioritized results
        
        Args:
            bed_path (str): Path to input BED file with methylation regions
        
        Returns:
            bool: True if processing succeeded, False otherwise
        """
        try:
            # 1. Map variants across regions
            mapped_variants = self.run_mapping(bed_path)
            
            if not mapped_variants:
                self.logger.error("No variants found during mapping")
                return False
            
            # 2. Score and prioritize variants
            meth_variants, unmeth_variants, rare_variants = self.score_variants(mapped_variants)  # Get rare_variants

            # 3. Write prioritized variant files
            self.write_prioritized_variants(meth_variants, unmeth_variants, rare_variants)  # Pass rare_variants
            
            # 4. Create overall statistics file
            self._write_overall_stats(mapped_variants)
            
            return True
        
        except Exception as e:
            self.logger.error(f"Comprehensive variant analysis failed: {str(e)}")
            return False

    def _write_overall_stats(self, mapped_variants: List[Dict]):
        """
        Write comprehensive statistics about the variant mapping and scoring
        
        Args:
            mapped_variants (List[Dict]): List of mapped variants
        """
        stats_file = f"{self.output_prefix}_analysis_stats.txt"
        
        try:
            with open(stats_file, 'w') as f:
                f.write("Variant Mapping and Scoring Analysis Statistics\n")
                f.write("==============================================\n\n")
                
                # Total variant counts
                f.write(f"Total Mapped Variants: {len(mapped_variants)}\n\n")
                
                # Variant type distribution
                type_counts = defaultdict(int)
                for variant in mapped_variants:
                    type_counts[variant['variant']['type']] += 1
                
                f.write("Variant Types Distribution:\n")
                for var_type, count in sorted(type_counts.items()):
                    f.write(f"  {var_type}: {count} ({count/len(mapped_variants)*100:.2f}%)\n")
                
                # Sample involvement statistics
                total_meth_samples = sum(v['total_meth'] for v in mapped_variants)
                total_unmeth_samples = sum(v['total_unmeth'] for v in mapped_variants)
                total_no_data_samples = sum(v['total_no_data'] for v in mapped_variants)
                
                f.write("\nSample Involvement:\n")
                f.write(f"  Methylated Samples: {total_meth_samples}\n")
                f.write(f"  Unmethylated Samples: {total_unmeth_samples}\n")
                f.write(f"  No Data Samples: {total_no_data_samples}\n")
            
            self.logger.info(f"Wrote analysis statistics to {stats_file}")
        
        except Exception as e:
            self.logger.error(f"Error writing analysis statistics: {str(e)}")

def main():
    """
    Main entry point for the Variant Methylation Mapper
    Handles command-line argument parsing and workflow execution
    """
    parser = argparse.ArgumentParser(
        description='Map and Score Variants Across Methylation Regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example Usage:
  python variant_methylation_mapper.py \\
    --bed methylation_regions.bed \\
    --haplotype H1 \\
    --output-prefix results/variant_analysis \\
    --small-vcf small_variants.vcf.gz \\
    --cnv-vcf cnvs.vcf.gz \\
    --sv-vcf structural_variants.vcf.gz \\
    --tr-vcf tandem_repeats.vcf.gz
        """
    )
    
    # Required arguments
    parser.add_argument('--bed', required=True,
                        help='Input BED file with methylation regions')
    parser.add_argument('--haplotype', required=True, 
                        choices=['H1', 'H2'],
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
        sys.exit(1)
    
    # Run complete workflow
    if not mapper.run(args.bed):
        mapper.logger.error("Variant mapping and scoring failed")
        sys.exit(1)
    
    mapper.logger.info("Variant Methylation Mapping and Scoring Completed Successfully")

if __name__ == '__main__':
    main()