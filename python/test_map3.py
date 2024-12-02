#!/usr/bin/env python3
"""
Optimized test script for mapping variants to H1/H2 methylation regions
Maps variants from VCF files and compares methylation patterns across haplotypes
"""
import argparse
import subprocess
import pandas as pd 
import os
import logging
from typing import Dict, List, Set, Tuple
from collections import defaultdict
import gzip
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
from tqdm import tqdm

class TestMapper:
    def __init__(self, output_prefix: str):
        self.output_prefix = output_prefix
        self.vcf_files = {}
        self.vcf_samples = {}
        self.sample_map = {}
        self.variant_counts = defaultdict(int)
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    def count_variants(self, vcf_file: str):
        with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                ref = fields[3]
                alt = fields[4]
                if len(ref) == 1 and len(alt) == 1:
                    self.variant_counts['SNP'] += 1
                elif len(ref) == 1 and len(alt) > 1:
                    self.variant_counts['Insertion'] += 1
                elif len(ref) > 1 and len(alt) == 1:
                    self.variant_counts['Deletion'] += 1
                else:
                    self.variant_counts['Complex'] += 1

    def log_variant_counts(self):
        logging.info('Variant counts by type:')
        for variant_type, count in self.variant_counts.items():
            logging.info(f'{variant_type}: {count}')

    def _setup_logging(self):
        """Configure logging"""
        os.makedirs(os.path.dirname(self.output_prefix), exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(f"{self.output_prefix}.log"),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def _normalize_sample(self, sample: str) -> str:
        """Normalize sample name with special handling for SPM276"""
        if sample in ('SPM276_1', 'SPM276-1'):
            return self.sample_map.get('SPM276_1', sample)
        return self.sample_map.get(sample.replace('-', '_'), sample)

    def _parse_samples(self, sample_str: str) -> Set[str]:
        """Parse sample string into normalized set"""
        if pd.isna(sample_str) or sample_str == '.':
            return set()
        return {self._normalize_sample(s.strip()) 
                for s in sample_str.replace(';', ',').split(',') if s.strip()}

    def load_vcfs(self, vcf_files: Dict[str, str]) -> bool:
        """Load and validate VCF files"""
        for vcf_type, path in vcf_files.items():
            if not path or not os.path.exists(path):
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
                                self.sample_map[s.replace('-', '_')] = s
                                self.sample_map[s] = s
                                if s in ('SPM276_1', 'SPM276-1'):
                                    self.sample_map['SPM276_1'] = s
                                    self.sample_map['SPM276-1'] = s
                            
                            self.logger.info(f"Loaded {vcf_type} VCF: {len(samples)} samples")
                            break
            except Exception as e:
                self.logger.error(f"Error loading {vcf_type} VCF: {str(e)}")
                
        return bool(self.vcf_files)

    def _get_variant_info(self, fields: List[str], vcf_type: str) -> Dict:
        """Extract basic variant information"""
        variant_id = fields[2]
        if variant_id == '.':
            variant_id = f"{vcf_type}_{fields[0]}_{fields[1]}"
            
        return {
            'variant_chr': fields[0],
            'variant_start': int(fields[1]),
            'variant_end': int(fields[1]) + len(fields[3]) - 1,
            'variant_id': variant_id,
            'type': vcf_type
        }

    def _process_genotypes(self, genotypes: Dict[str, str], hap_idx: int,
                          meth_samples: Set[str], unmeth_samples: Set[str]) -> Tuple[List[str], List[str]]:
        """Process genotypes for a single haplotype"""
        meth_vars, unmeth_vars = [], []
        
        for sample, gt_str in genotypes.items():
            if gt_str in {'./.', '.|.'}:
                continue
                
            gt = gt_str.split(':')[0]
            if '|' not in gt:
                continue
                
            alleles = gt.split('|')
            if len(alleles) != 2 or int(alleles[hap_idx]) == 0:
                continue
                
            if sample in meth_samples:
                meth_vars.append(f"{sample}:{alleles[hap_idx]}")
            elif sample in unmeth_samples:
                unmeth_vars.append(f"{sample}:{alleles[hap_idx]}")
                
        return meth_vars, unmeth_vars

    def process_region(self, h1_region: pd.Series, h2_region: pd.Series) -> List[Dict]:
        """Process variants in a region pair"""
        results = []
        
        # Get sample sets
        h1_meth = self._parse_samples(h1_region.methylated_samples)
        h1_unmeth = self._parse_samples(h1_region.unmethylated_samples) 
        h2_meth = self._parse_samples(h2_region.methylated_samples)
        h2_unmeth = self._parse_samples(h2_region.unmethylated_samples)

        for vcf_type, vcf_path in self.vcf_files.items():
            try:
                # Get variants in region
                cmd = ['tabix', vcf_path, f"{h1_region.chrom}:{h1_region.start}-{h1_region.end}"]
                output = subprocess.check_output(cmd, text=True)

                for line in output.splitlines():
                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        continue

                    # Get variant info and genotypes
                    variant = self._get_variant_info(fields, vcf_type)
                    genotypes = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
                    
                    # Process each haplotype
                    h1_meth_vars, h1_unmeth_vars = self._process_genotypes(
                        genotypes, 0, h1_meth, h1_unmeth)
                    h2_meth_vars, h2_unmeth_vars = self._process_genotypes(
                        genotypes, 1, h2_meth, h2_unmeth)

                    if any([h1_meth_vars, h1_unmeth_vars, h2_meth_vars, h2_unmeth_vars]):
                        variant.update({
                            'chrom': h1_region.chrom,
                            'start': h1_region.start,
                            'end': h1_region.end,
                            'h1_meth_frac': f"{len(h1_meth_vars)}/{len(h1_meth)}",
                            'h1_unmeth_frac': f"{len(h1_unmeth_vars)}/{len(h1_unmeth)}",
                            'h2_meth_frac': f"{len(h2_meth_vars)}/{len(h2_meth)}",
                            'h2_unmeth_frac': f"{len(h2_unmeth_vars)}/{len(h2_unmeth)}",
                            'h1_meth_samples': ','.join(h1_meth_vars) or '.',
                            'h1_unmeth_samples': ','.join(h1_unmeth_vars) or '.',
                            'h2_meth_samples': ','.join(h2_meth_vars) or '.',
                            'h2_unmeth_samples': ','.join(h2_unmeth_vars) or '.'
                        })
                        results.append(variant)

            except subprocess.CalledProcessError:
                pass
            except Exception as e:
                self.logger.error(f"Error processing variants: {str(e)}")

        return results

    def test(self, h1_bed: str, h2_bed: str, sample_size: int = 50) -> bool:
        """Run test mapping on sample regions"""
        try:
            # Load regions
            h1_df = pd.read_csv(h1_bed, sep='\t')
            h2_df = pd.read_csv(h2_bed, sep='\t')
            
            # Find matching regions
            h1_df['coord'] = h1_df['chrom'] + '_' + h1_df['start'].astype(str)
            h2_df['coord'] = h2_df['chrom'] + '_' + h2_df['start'].astype(str)
            shared = set(h1_df['coord']) & set(h2_df['coord'])
            
            if not shared:
                self.logger.error("No matching regions found between H1 and H2")
                return False
                
            # Sample regions
            test_coords = pd.Series(list(shared)).sample(n=min(sample_size, len(shared)))
            h1_test = h1_df[h1_df['coord'].isin(test_coords)]
            h2_test = h2_df[h2_df['coord'].isin(test_coords)]
            
            # Process regions
            results = []
            with ThreadPoolExecutor(max_workers=8) as executor:
                futures = [
                    executor.submit(self.process_region, h1_row, h2_row)
                    for h1_row, h2_row in zip(h1_test.itertuples(), h2_test.itertuples())
                ]
                
                for future in tqdm(futures, desc="Processing regions"):
                    try:
                        results.extend(future.result())
                    except Exception as e:
                        self.logger.error(f"Error in region: {str(e)}")

            if not results:
                self.logger.error("No variants found")
                return False

            # Write results
            pd.DataFrame(results).to_csv(self.output_prefix, sep='\t', index=False)
            self.logger.info(f"Found {len(results)} variants in {len(test_coords)} regions")
            return True

        except Exception as e:
            self.logger.error(f"Test failed: {str(e)}")
            return False

def main():
    parser = argparse.ArgumentParser(
        description='Test variant mapping across H1/H2 regions'
    )
    parser.add_argument('--h1-bed', required=True,
                       help='H1 methylation regions BED file')
    parser.add_argument('--h2-bed', required=True,
                       help='H2 methylation regions BED file')
    parser.add_argument('--output', required=True,
                       help='Output file path')
    parser.add_argument('--small-vcf',
                       help='Small variants VCF')
    parser.add_argument('--cnv-vcf',
                       help='CNV variants VCF')
    parser.add_argument('--sv-vcf',
                       help='Structural variants VCF')
    parser.add_argument('--sample-size', type=int, default=50,
                       help='Number of random regions to test')

    args = parser.parse_args()

    mapper = TestMapper(args.output)
    if not mapper.load_vcfs({
        'small': args.small_vcf,
        'cnv': args.cnv_vcf,
        'sv': args.sv_vcf
    }):
        mapper.logger.error("No valid VCF files provided")
        exit(1)

    if not mapper.test(args.h1_bed, args.h2_bed, args.sample_size):
        mapper.logger.error("Test failed")
        exit(1)

if __name__ == '__main__':
    main()