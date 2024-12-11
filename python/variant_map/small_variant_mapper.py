#!/usr/bin/env python3
"""
simple_variant_mapper.py - Maps variants from VCF to methylation regions
Uses basic file handling and tabix for VCF parsing.
"""

import argparse
import subprocess
import pandas as pd
import os
import logging
from typing import Dict, List, Set, Optional, Iterator
import gzip
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp

class VariantMapper:
    def __init__(self, output_file: str):
        self.output_file = output_file
        self.vcf_path = None
        self.sample_name_map = {}
        self.vcf_samples = []
        self.setup_logging()
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

    def setup_logging(self):
        """Configure logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)

    def _parse_vcf_header(self, vcf_path: str) -> bool:
        """Parse VCF header to get samples"""
        try:
            with gzip.open(vcf_path, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        if line.startswith('#CHROM'):
                            # Get samples from header line
                            fields = line.strip().split('\t')
                            self.vcf_samples = fields[9:]  # Samples start at column 9
                            self._build_sample_map(self.vcf_samples)
                            return True
                    else:
                        break
            return False
        except Exception as e:
            self.logger.error(f"Error parsing VCF header: {str(e)}")
            return False

    def _build_sample_map(self, samples: List[str]):
        """Build sample name mapping with SPM276 handling"""
        for sample in samples:
            # Handle SPM276 case specifically
            if sample in ('SPM276_1', 'SPM276-1'):
                self.sample_name_map['SPM276_1'] = sample
                self.sample_name_map['SPM276-1'] = sample
            else:
                # Handle other sample names
                self.sample_name_map[sample] = sample
                if '-' in sample:
                    self.sample_name_map[sample.replace('-', '_')] = sample
                if '_' in sample:
                    self.sample_name_map[sample.replace('_', '-')] = sample

    def load_vcf(self, vcf_path: str) -> bool:
        """Load and validate VCF file"""
        if not vcf_path or not os.path.exists(vcf_path):
            self.logger.error(f"VCF file not found: {vcf_path}")
            return False

        try:
            # Check if index exists
            if not os.path.exists(vcf_path + '.tbi'):
                self.logger.error(f"VCF index not found: {vcf_path}.tbi")
                return False

            # Parse header
            if not self._parse_vcf_header(vcf_path):
                self.logger.error("Failed to parse VCF header")
                return False

            self.vcf_path = vcf_path
            self.logger.info(f"Loaded VCF file: {vcf_path}")
            self.logger.info(f"Found {len(self.vcf_samples)} samples")
            return True

        except Exception as e:
            self.logger.error(f"Error loading VCF file: {str(e)}")
            return False

    def _normalize_sample_name(self, sample: str) -> str:
        """Get normalized sample name from mapping"""
        if sample in ('SPM276_1', 'SPM276-1'):
            return self.sample_name_map.get('SPM276_1', sample)
        return self.sample_name_map.get(sample, sample)

    def _parse_samples(self, sample_str: str) -> Set[str]:
        """Parse sample string handling both comma and semicolon separators"""
        if pd.isna(sample_str) or sample_str == '.':
            return set()
        
        if ',' in sample_str:
            samples = {s.strip() for s in sample_str.split(',') if s.strip()}
        elif ';' in sample_str:
            samples = {s.strip() for s in sample_str.split(';') if s.strip()}
        else:
            samples = {sample_str.strip()}
            
        return {self._normalize_sample_name(s) for s in samples}

    def _fetch_variants(self, chrom: str, start: int, end: int) -> Iterator[Dict]:
        """Fetch variants from VCF using tabix"""
        try:
            # Use tabix to fetch region
            cmd = ['tabix', self.vcf_path, f'{chrom}:{start}-{end}']
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            
            if stderr:
                self.logger.error(f"Tabix error: {stderr.decode()}")
                return

            # Parse variants
            for line in stdout.decode().split('\n'):
                if not line.strip():
                    continue
                    
                fields = line.split('\t')
                if len(fields) < 10:  # Must have at least mandatory fields + 1 sample
                    continue

                # Parse basic variant info
                var_info = {
                    'chrom': fields[0],
                    'pos': int(fields[1]),
                    'id': fields[2],
                    'ref': fields[3],
                    'alt': fields[4],
                    'genotypes': dict(zip(self.vcf_samples, fields[9:]))
                }
                yield var_info

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Tabix command failed: {e}")
        except Exception as e:
            self.logger.error(f"Error fetching variants: {str(e)}")

    def _parse_genotype(self, gt_str: str) -> Optional[str]:
        """Parse genotype string to get basic GT field"""
        if not gt_str or gt_str == './.':
            return None
        try:
            gt = gt_str.split(':')[0]  # Get GT field
            if any(int(a) > 0 for a in gt.replace('|', '/').split('/') if a.isdigit()):
                return gt
        except:
            pass
        return None

    def process_region(self, region: pd.Series) -> List[Dict]:
        """Process variants in a region"""
        variants = []
        
        try:
            # Get methylated and unmethylated samples
            meth_samples = self._parse_samples(region['methylated_samples'])
            unmeth_samples = self._parse_samples(region['unmethylated_samples'])

            # Process variants in region
            for variant in self._fetch_variants(
                region['chrom'], 
                int(region['start']), 
                int(region['end'])
            ):
                meth_carriers = []
                unmeth_carriers = []

                # Check samples
                for sample, gt_str in variant['genotypes'].items():
                    gt = self._parse_genotype(gt_str)
                    if gt:
                        sample = self._normalize_sample_name(sample)
                        if sample in meth_samples:
                            meth_carriers.append(f"{sample}:{gt}")
                        elif sample in unmeth_samples:
                            unmeth_carriers.append(f"{sample}:{gt}")

                # Record variant if found in relevant samples
                if meth_carriers or unmeth_carriers:
                    variants.append({
                        'chrom': region['chrom'],
                        'start': region['start'],
                        'end': region['end'],
                        'variant_id': variant['id'],
                        'position': variant['pos'],
                        'ref': variant['ref'],
                        'alt': variant['alt'],
                        'meth_samples': ','.join(meth_carriers) if meth_carriers else '.',
                        'unmeth_samples': ','.join(unmeth_carriers) if unmeth_carriers else '.'
                    })

        except Exception as e:
            self.logger.error(f"Error processing region {region['chrom']}:{region['start']}-{region['end']}: {str(e)}")
            
        return variants

    def process_bed(self, bed_path: str) -> Optional[List[Dict]]:
        """Process BED file"""
        self.logger.info(f"Processing BED file: {bed_path}")
        
        try:
            results = []
            df = pd.read_csv(bed_path, sep='\t')
            
            # Process in batches using ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=4) as executor:
                futures = []
                
                for _, region in df.iterrows():
                    futures.append(executor.submit(self.process_region, region))
                
                for future in futures:
                    try:
                        variants = future.result()
                        results.extend(variants)
                    except Exception as e:
                        self.logger.error(f"Error processing region: {str(e)}")
            
            self.logger.info(f"Found {len(results)} variants")
            return results
            
        except Exception as e:
            self.logger.error(f"Error processing BED file: {str(e)}")
            return None

    def write_results(self, results: List[Dict]) -> bool:
        """Write results to output file"""
        if not results:
            self.logger.warning("No results to write")
            return False
            
        try:
            headers = ['chrom', 'start', 'end', 'variant_id', 'position',
                      'ref', 'alt', 'meth_samples', 'unmeth_samples']
            
            with open(self.output_file, 'w') as f:
                f.write('\t'.join(headers) + '\n')
                for result in results:
                    row = [str(result.get(h, '.')) for h in headers]
                    f.write('\t'.join(row) + '\n')
                    
            self.logger.info(f"Wrote {len(results)} variants to {self.output_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing results: {str(e)}")
            return False

def main():
    parser = argparse.ArgumentParser(description='Map variants to methylation regions')
    parser.add_argument('--bed', required=True, help='Input BED file')
    parser.add_argument('--vcf', required=True, help='Input VCF file (must be bgzipped and indexed)')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()
    
    mapper = VariantMapper(args.output)
    
    if not mapper.load_vcf(args.vcf):
        exit(1)
    
    results = mapper.process_bed(args.bed)
    if results:
        if not mapper.write_results(results):
            exit(1)
    else:
        mapper.logger.error("No results generated")
        exit(1)

if __name__ == '__main__':
    main()