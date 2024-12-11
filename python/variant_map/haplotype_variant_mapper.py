#!/usr/bin/env python3
"""
variant_mapper.py - Production-ready variant mapper for methylation regions
Maps variants from VCF files (SNPs, CNVs, SVs, TRs) to methylation regions with haplotype specificity
"""
import argparse
import subprocess
import pandas as pd
import os
import logging
from typing import Dict, List, Set, Optional, Iterator
import gzip
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from datetime import datetime
from functools import lru_cache
import numpy as np
from tqdm import tqdm

class VariantMapper:
    def __init__(self, output_file: str, haplotype: str):
        self.output_file = output_file
        self.haplotype = haplotype
        self.vcf_files = {}
        self.vcf_samples = {}
        self.sample_map = {}
        self.region_cache = {}
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Setup logging
        log_file = os.path.join(os.path.dirname(output_file), 
                               f'mapper_{haplotype}_{datetime.now():%Y%m%d_%H%M%S}.log')
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
            'small': {
                'phased': True,
                'format_fields': ['GT', 'AD', 'DP', 'GQ'],
                'description': 'DeepVariant SNPs/indels'
            },
            'cnv': {
                'phased': False,
                'format_fields': ['GT', 'CN', 'CNQ'],
                'description': 'HiFiCNV variants'
            },
            'sv': {
                'phased': True,
                'format_fields': ['GT', 'DR', 'DV'],
                'description': 'PBSV structural variants'
            },
            'tr': {
                'phased': False,
                'format_fields': ['GT', 'RU', 'RC'],
                'description': 'TRGT tandem repeats'
            }
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
                                norm_s = s.replace('-','_')
                                self.sample_map[norm_s] = s
                                self.sample_map[s] = s
                                # Special case for SPM276
                                if s in ('SPM276_1', 'SPM276-1'):
                                    self.sample_map['SPM276_1'] = s
                                    self.sample_map['SPM276-1'] = s
                                    
                            success = True
                            self.logger.info(f"Loaded {vcf_type} VCF: {len(samples)} samples")
                            break
                            
            except Exception as e:
                self.logger.error(f"Error loading {vcf_type} VCF: {str(e)}")
                
        return success

    @lru_cache(maxsize=10000)
    def _get_samples(self, sample_str: str) -> Set[str]:
        """Parse and normalize sample names with caching"""
        if not sample_str or sample_str in {'.', 'nan'}:
            return set()
        return {self.sample_map.get(s.strip(), s.strip()) 
                for s in sample_str.replace(';',',').split(',') if s.strip()}

    def _parse_genotype(self, gt_str: str, vcf_type: str) -> Optional[str]:
        """Parse genotype string based on variant type"""
        if not gt_str or gt_str in {'./.', '.|.'}:
            return None
            
        fields = gt_str.split(':')
        gt = fields[0]
        
        try:
            if vcf_type == 'cnv':
                cn = next((f for f in fields if f.startswith('CN=')), None)
                if cn:
                    cn_val = int(cn.split('=')[1])
                    if cn_val != 2:  # Non-reference CN
                        return f"CN={cn_val}"
                        
            elif vcf_type == 'tr':
                ru = next((f for f in fields if f.startswith('RU=')), None)
                rc = next((f for f in fields if f.startswith('RC=')), None)
                if ru and rc:
                    return f"{ru.split('=')[1]}[{rc.split('=')[1]}]"
                    
            elif self.variant_configs[vcf_type]['phased']:
                if '|' in gt:
                    hap_idx = 0 if self.haplotype == 'H1' else 1
                    alleles = gt.split('|')
                    if len(alleles) == 2 and int(alleles[hap_idx]) > 0:
                        # Add depth info if available
                        dp = next((f for f in fields if f.startswith('DP=')), None)
                        ad = next((f for f in fields if f.startswith('AD=')), None)
                        if dp and ad:
                            return f"{gt}:{ad}:{dp}"
                        return gt
                        
            else:  # Unphased variants
                if any(int(a) > 0 for a in gt.replace('|','/').split('/') 
                      if a.isdigit()):
                    return gt
                    
        except Exception:
            pass
            
        return None

    def _format_alt(self, variant: Dict, vcf_type: str) -> str:
        """Format alternative allele based on variant type"""
        if vcf_type == 'small':
            return variant['alt']
            
        if vcf_type == 'sv':
            sv_type = variant['info'].get('SVTYPE', 'SV')
            sv_len = variant['info'].get('SVLEN', '')
            sv_len = f":{abs(int(sv_len))}" if sv_len else ''
            return f"<{sv_type}{sv_len}>"
            
        if vcf_type == 'cnv':
            cn = variant['info'].get('CN', '?')
            return f"<CNV:{cn}>"
            
        if vcf_type == 'tr':
            ru = variant['info'].get('RU', variant['ref'])
            rc = variant['info'].get('RC', '?')
            return f"{ru}[{rc}]"
            
        return variant['alt']

    def process_region(self, region: pd.Series) -> List[Dict]:
        """Process variants in a genomic region"""
        cache_key = f"{region['chrom']}:{region['start']}-{region['end']}"
        if cache_key in self.region_cache:
            return self.region_cache[cache_key]
            
        variants = []
        meth = self._get_samples(region['methylated_samples'])
        unmeth = self._get_samples(region['unmethylated_samples'])
        
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
                        
                    # Parse INFO field
                    info = dict(item.split('=') if '=' in item else (item, True)
                              for item in fields[7].split(';')) if fields[7] != '.' else {}
                              
                    # Get genotypes
                    gts = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
                    meth_vars = []
                    unmeth_vars = []
                    
                    # Process methylated and unmethylated samples
                    for samples, carriers in [(meth, meth_vars), (unmeth, unmeth_vars)]:
                        for sample in samples:
                            if sample not in gts:
                                continue
                                
                            gt = self._parse_genotype(gts[sample], vcf_type)
                            if gt:
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
                            'alt': self._format_alt({'alt': fields[4], 'info': info}, vcf_type),
                            'num_meth': len(meth_vars),
                            'num_unmeth': len(unmeth_vars),
                            'meth_samples': ','.join(meth_vars) if meth_vars else '.',
                            'unmeth_samples': ','.join(unmeth_vars) if unmeth_vars else '.'
                        }
                        variants.append(variant)
                        
            except subprocess.CalledProcessError:
                continue
            except Exception as e:
                self.logger.error(f"Error processing {vcf_type} variants: {str(e)}")
                
        self.region_cache[cache_key] = variants
        return variants

    def process_bed(self, bed_path: str) -> List[Dict]:
        """Process BED file in parallel with progress tracking"""
        self.logger.info(f"Reading BED file: {bed_path}")
        
        results = []
        chunk_size = 1000
        max_workers = min(mp.cpu_count() * 2, 16)
        
        # Process in chunks with progress bar
        chunks = pd.read_csv(bed_path, sep='\t', chunksize=chunk_size)
        total_regions = sum(1 for _ in pd.read_csv(bed_path, sep='\t', chunksize=chunk_size))
        
        with tqdm(total=total_regions, desc="Processing regions") as pbar:
            for chunk in chunks:
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    futures = [executor.submit(self.process_region, region) 
                             for _, region in chunk.iterrows()]
                    
                    for f in futures:
                        try:
                            variants = f.result()
                            results.extend(variants)
                            pbar.update(1)
                        except Exception as e:
                            self.logger.error(f"Error processing region: {str(e)}")
                            pbar.update(1)
                            
        self.logger.info(f"Found {len(results)} total variants")
        return results

    def write_results(self, results: List[Dict]) -> bool:
        """Write results with detailed statistics"""
        if not results:
            return False
            
        # Write main results
        headers = [
            'chrom', 'start', 'end', 'variant_id', 'type', 'ref', 'alt',
            'num_meth', 'num_unmeth', 'meth_samples', 'unmeth_samples'
        ]
        
        try:
            with open(self.output_file, 'w', buffering=16*1024*1024) as f:
                f.write('\t'.join(headers) + '\n')
                for r in results:
                    f.write('\t'.join(str(r.get(h, '.')) for h in headers) + '\n')
                    
            # Generate statistics
            stats_file = f"{os.path.splitext(self.output_file)[0]}_stats.txt"
            with open(stats_file, 'w') as f:
                f.write(f"Variant Mapping Statistics ({self.haplotype})\n")
                f.write("=" * 50 + "\n\n")
                
                # Variant type counts
                type_counts = {}
                for r in results:
                    type_counts[r['type']] = type_counts.get(r['type'], 0) + 1
                
                f.write("Variant Counts by Type:\n")
                for vtype, count in type_counts.items():
                    f.write(f"{vtype:10} : {count:8,d}\n")
                f.write(f"{'Total':10} : {len(results):8,d}\n\n")
                
                # Sample counts
                all_samples = set()
                meth_samples = set()
                unmeth_samples = set()
                
                for r in results:
                    if r['meth_samples'] != '.':
                        curr_samples = {s.split(':')[0] for s in r['meth_samples'].split(',')}
                        all_samples.update(curr_samples)
                        meth_samples.update(curr_samples)
                    if r['unmeth_samples'] != '.':
                        curr_samples = {s.split(':')[0] for s in r['unmeth_samples'].split(',')}
                        all_samples.update(curr_samples)
                        unmeth_samples.update(curr_samples)
                
                f.write("Sample Statistics:\n")
                f.write(f"Total unique samples      : {len(all_samples):,d}\n")
                f.write(f"Samples with methylation  : {len(meth_samples):,d}\n")
                f.write(f"Samples with unmethylation: {len(unmeth_samples):,d}\n")
                f.write(f"Samples with both states  : {len(meth_samples & unmeth_samples):,d}\n")
                
            self.logger.info(f"Results written to: {self.output_file}")
            self.logger.info(f"Statistics written to: {stats_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing results: {str(e)}")
            return False

def main():
    parser = argparse.ArgumentParser(description='Map variants to methylation regions')
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
    parser.add_argument('--threads', type=int,
                       help='Number of processing threads (default: auto)')
    
    args = parser.parse_args()
    
    # Initialize mapper
    mapper = VariantMapper(args.output, args.haplotype)
    
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
    
    # Process variants
    try:
        results = mapper.process_bed(args.bed)
        if not results or not mapper.write_results(results):
            mapper.logger.error("Error processing variants")
            exit(1)
    except Exception as e:
        mapper.logger.error(f"Processing failed: {str(e)}")
        exit(1)

if __name__ == '__main__':
    main()