#!/usr/bin/env python3
"""
quick_test_mapper.py - Fast variant mapping test using 50 random regions
Tests variant mapping for small variants, CNVs, SVs, and tandem repeats against methylation data
"""
import argparse
import subprocess
import pandas as pd
import os
import logging
from typing import Dict, List, Set
import gzip
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from datetime import datetime

class QuickTestMapper:
    def __init__(self, output_file: str, haplotype: str):
        self.output_file = output_file
        self.haplotype = haplotype
        self.vcf_files = {}
        self.vcf_samples = {}
        self.sample_map = {}
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Setup logging
        log_file = os.path.join(os.path.dirname(output_file), 
                               f'mapper_test_{datetime.now():%Y%m%d_%H%M%S}.log')
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # Variant configurations
        self.variant_configs = {
            'small': {'phased': True, 'description': 'DeepVariant SNPs/indels'},
            'cnv': {'phased': False, 'description': 'HiFiCNV variants'},
            'sv': {'phased': True, 'description': 'PBSV structural variants'},
            'tr': {'phased': False, 'description': 'TRGT tandem repeats'}
        }

    def load_vcfs(self, vcf_paths: Dict[str, str]) -> bool:
        """Load and validate VCF files"""
        success = False
        for vcf_type, path in vcf_paths.items():
            if not path or not os.path.exists(path):
                self.logger.warning(f"Skipping {vcf_type}: file not found")
                continue
                
            if not os.path.exists(f"{path}.tbi"):
                self.logger.warning(f"Skipping {vcf_type}: index not found")
                continue
                
            try:
                with gzip.open(path, 'rt') as f:
                    for line in f:
                        if line.startswith('#CHROM'):
                            samples = line.strip().split('\t')[9:]
                            self.vcf_samples[vcf_type] = samples
                            self.vcf_files[vcf_type] = path
                            # Build sample name mapping
                            for s in samples:
                                self.sample_map[s.replace('-','_')] = s
                                self.sample_map[s] = s
                                # Special case for SPM276
                                if s in ('SPM276_1', 'SPM276-1'):
                                    self.sample_map['SPM276_1'] = s
                                    self.sample_map['SPM276-1'] = s
                            success = True
                            self.logger.info(f"Loaded {vcf_type} VCF with {len(samples)} samples")
                            break
            except Exception as e:
                self.logger.error(f"Error loading {vcf_type} VCF: {str(e)}")
                
        return success

    def process_region(self, region: pd.Series) -> List[Dict]:
        """Process variants in a single region"""
        variants = []
        
        # Get methylated and unmethylated samples
        try:
            meth = set(s.strip() for s in str(region['methylated_samples']).split(',') 
                      if s.strip() not in ['.', 'nan'])
            unmeth = set(s.strip() for s in str(region['unmethylated_samples']).split(',') 
                        if s.strip() not in ['.', 'nan'])
        except Exception as e:
            self.logger.error(f"Error parsing samples in region {region['chrom']}:{region['start']}: {e}")
            return []

        # Process each VCF file
        for vcf_type, vcf_path in self.vcf_files.items():
            try:
                # Fetch variants using tabix
                cmd = ['tabix', vcf_path, f"{region['chrom']}:{region['start']}-{region['end']}"]
                output = subprocess.check_output(cmd, text=True)
                
                # Process each variant
                for line in output.splitlines():
                    fields = line.split('\t')
                    if len(fields) < 10:
                        continue
                        
                    # Get sample genotypes
                    gts = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
                    meth_vars = []
                    unmeth_vars = []
                    
                    # Check methylated and unmethylated samples
                    for samples, carriers in [(meth, meth_vars), (unmeth, unmeth_vars)]:
                        for sample in samples:
                            if sample not in gts:
                                continue
                                
                            gt = gts[sample].split(':')[0]
                            if gt in ['./.', '.|.']:
                                continue
                                
                            # Handle phased/unphased variants
                            if self.variant_configs[vcf_type]['phased']:
                                if '|' in gt:
                                    hap_idx = 0 if self.haplotype == 'H1' else 1
                                    alleles = gt.split('|')
                                    if len(alleles) == 2 and int(alleles[hap_idx]) > 0:
                                        carriers.append(f"{sample}:{gt}")
                            else:
                                if any(int(a) > 0 for a in gt.replace('|','/').split('/') 
                                      if a.isdigit()):
                                    carriers.append(f"{sample}:{gt}")
                    
                    # Record variant if it appears in any samples
                    if meth_vars or unmeth_vars:
                        variant = {
                            'chrom': fields[0],
                            'start': region['start'],
                            'end': region['end'],
                            'variant_id': fields[2] if fields[2] != '.' else f"{vcf_type}_{fields[0]}_{fields[1]}",
                            'type': vcf_type,
                            'ref': fields[3],
                            'alt': fields[4],
                            'num_meth': len(meth_vars),
                            'num_unmeth': len(unmeth_vars),
                            'meth_samples': ','.join(meth_vars) if meth_vars else '.',
                            'unmeth_samples': ','.join(unmeth_vars) if unmeth_vars else '.'
                        }
                        variants.append(variant)
                        
            except subprocess.CalledProcessError:
                self.logger.debug(f"No variants found in {vcf_type} for region {region['chrom']}:{region['start']}")
            except Exception as e:
                self.logger.error(f"Error processing {vcf_type} variants: {str(e)}")
                
        return variants

    def run_test(self, bed_path: str) -> List[Dict]:
        """Run test on sample of regions"""
        self.logger.info(f"Reading BED file: {bed_path}")
        df = pd.read_csv(bed_path, sep='\t')
        self.logger.info(f"Found {len(df)} total regions")
        
        # Sample 50 random regions
        sample_size = min(50, len(df))
        sample_regions = df.sample(n=sample_size)
        self.logger.info(f"Randomly selected {sample_size} regions for testing")
        
        results = []
        max_workers = min(mp.cpu_count(), 8)
        
        # Process regions in parallel
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self.process_region, region) 
                      for _, region in sample_regions.iterrows()]
            
            for f in futures:
                try:
                    variants = f.result()
                    results.extend(variants)
                except Exception as e:
                    self.logger.error(f"Error processing region: {str(e)}")
                    
        self.logger.info(f"Found {len(results)} variants across {sample_size} regions")
        return results

    def write_results(self, results: List[Dict]) -> bool:
        """Write results to TSV file"""
        if not results:
            self.logger.warning("No results to write")
            return False
            
        headers = [
            'chrom', 'start', 'end', 'variant_id', 'type', 'ref', 'alt',
            'num_meth', 'num_unmeth', 'meth_samples', 'unmeth_samples'
        ]
            
        try:
            with open(self.output_file, 'w') as f:
                f.write('\t'.join(headers) + '\n')
                for r in results:
                    f.write('\t'.join(str(r.get(h, '.')) for h in headers) + '\n')
                    
            self.logger.info(f"Results written to: {self.output_file}")
            
            # Log variant type counts
            type_counts = {}
            for r in results:
                type_counts[r['type']] = type_counts.get(r['type'], 0) + 1
            self.logger.info("Variant counts by type:")
            for vtype, count in type_counts.items():
                self.logger.info(f"  {vtype}: {count}")
                
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing results: {str(e)}")
            return False

def main():
    parser = argparse.ArgumentParser(
        description='Quick test variant mapping on sample of 50 regions'
    )
    parser.add_argument('--bed', required=True,
                       help='Input BED file with methylation regions')
    parser.add_argument('--haplotype', required=True, choices=['H1', 'H2'],
                       help='Haplotype to process')
    parser.add_argument('--output', required=True,
                       help='Output file path')
    parser.add_argument('--small-vcf',
                       help='Small variants VCF (DeepVariant)')
    parser.add_argument('--cnv-vcf',
                       help='Copy number variants VCF (HiFiCNV)')
    parser.add_argument('--sv-vcf',
                       help='Structural variants VCF (PBSV)')
    parser.add_argument('--tr-vcf',
                       help='Tandem repeats VCF (TRGT)')
    
    args = parser.parse_args()
    
    # Initialize mapper
    mapper = QuickTestMapper(args.output, args.haplotype)
    
    # Load VCF files
    vcf_files = {
        'small': args.small_vcf,
        'cnv': args.cnv_vcf,
        'sv': args.sv_vcf,
        'tr': args.tr_vcf
    }
    
    if not mapper.load_vcfs(vcf_files):
        mapper.logger.error("No valid VCF files provided")
        exit(1)
    
    # Run test
    try:
        results = mapper.run_test(args.bed)
        if not results or not mapper.write_results(results):
            mapper.logger.error("Error processing variants")
            exit(1)
    except Exception as e:
        mapper.logger.error(f"Test failed: {str(e)}")
        exit(1)

if __name__ == '__main__':
    main()