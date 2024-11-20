#!/usr/bin/env python3
"""
simplified_variant_mapper.py - Maps variants from VCF files to methylation regions
Supports mapping of small variants, CNVs, SVs, and tandem repeats.
Includes special handling for SPM276 sample naming conventions.
"""

import argparse
from pysam import VariantFile
import pandas as pd
import os
import logging
from typing import Dict, List, Set, Optional
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp

class VariantMapper:
    def __init__(self, output_file: str):
        """Initialize the variant mapper"""
        self.output_file = output_file
        self.vcf_handlers = {}
        self.sample_name_map = {}
        self.setup_logging()
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

    def setup_logging(self):
        """Configure logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)

    def _build_sample_map(self, vcf_samples: Set[str]):
        """Build sample name mapping with SPM276 handling"""
        for sample in vcf_samples:
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

    def add_variant_file(self, file_type: str, file_path: str) -> bool:
        """Add a VCF file for processing"""
        if not file_path or not os.path.exists(file_path):
            self.logger.error(f"{file_type} file not found: {file_path}")
            return False
            
        try:
            vcf = VariantFile(file_path)
            self._build_sample_map(set(vcf.header.samples))
            self.vcf_handlers[file_type] = vcf
            self.logger.info(f"Loaded {file_type} file: {file_path}")
            return True
        except Exception as e:
            self.logger.error(f"Error loading {file_type} file: {str(e)}")
            return False

    def _normalize_sample_name(self, sample: str) -> str:
        """Get normalized sample name from mapping"""
        # Handle SPM276 specifically
        if sample in ('SPM276_1', 'SPM276-1'):
            return self.sample_name_map.get('SPM276_1', sample)
        
        # Try direct lookup first
        if sample in self.sample_name_map:
            return self.sample_name_map[sample]
        
        # Try variant with swapped separator
        variant = sample.replace('-', '_') if '-' in sample else sample.replace('_', '-')
        return self.sample_name_map.get(variant, sample)

    def _parse_samples(self, sample_str: str) -> Set[str]:
        """Parse sample string handling both comma and semicolon separators"""
        if pd.isna(sample_str) or sample_str == '.':
            return set()
        
        # Try comma first, then semicolon
        if ',' in sample_str:
            samples = {s.strip() for s in sample_str.split(',') if s.strip()}
        elif ';' in sample_str:
            samples = {s.strip() for s in sample_str.split(';') if s.strip()}
        else:
            samples = {sample_str.strip()}
            
        return {self._normalize_sample_name(s) for s in samples}

    def process_region(self, region: pd.Series, vcf_type: str, vcf: VariantFile) -> List[Dict]:
        """Process variants in a region"""
        variants = []
        
        # Get methylated and unmethylated samples with normalization
        meth_samples = self._parse_samples(region['methylated_samples'])
        unmeth_samples = self._parse_samples(region['unmethylated_samples'])

        try:
            # Process each variant in the region
            for variant in vcf.fetch(region['chrom'], region['start'], region['end']):
                # Track which samples have the variant
                meth_carriers = []
                unmeth_carriers = []
                
                # Check methylated samples
                for sample in meth_samples:
                    if sample in variant.samples:
                        gt = variant.samples[sample]['GT']
                        if gt and any(g is not None and g > 0 for g in gt):
                            # Convert back to original format for output
                            orig_sample = sample.replace('_', '-')
                            meth_carriers.append(f"{orig_sample}:{gt[0]}/{gt[1]}")

                # Check unmethylated samples
                for sample in unmeth_samples:
                    if sample in variant.samples:
                        gt = variant.samples[sample]['GT']
                        if gt and any(g is not None and g > 0 for g in gt):
                            # Convert back to original format for output
                            orig_sample = sample.replace('_', '-')
                            unmeth_carriers.append(f"{orig_sample}:{gt[0]}/{gt[1]}")
                
                # Only record variant if it appears in our samples of interest
                if meth_carriers or unmeth_carriers:
                    variants.append({
                        'chrom': region['chrom'],
                        'start': region['start'],
                        'end': region['end'],
                        'variant_id': variant.id or f"{variant.chrom}:{variant.pos}",
                        'position': variant.pos,
                        'type': vcf_type,
                        'ref': variant.ref,
                        'alt': ','.join(str(a) for a in variant.alts) if variant.alts else '.',
                        'meth_samples': ','.join(meth_carriers) if meth_carriers else '.',
                        'unmeth_samples': ','.join(unmeth_carriers) if unmeth_carriers else '.'
                    })
                    
        except Exception as e:
            self.logger.error(f"Error processing region {region['chrom']}:{region['start']}-{region['end']}: {str(e)}")
            
        return variants

    def process_bed(self, bed_path: str) -> Optional[List[Dict]]:
        """Process BED file with parallel processing"""
        self.logger.info(f"Processing BED file: {bed_path}")
        
        try:
            results = []
            df = pd.read_csv(bed_path, sep='\t')
            
            # Use ThreadPoolExecutor for parallel processing
            max_workers = min(mp.cpu_count(), 8)
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = []
                
                # Submit tasks
                for _, region in df.iterrows():
                    for vcf_type, vcf in self.vcf_handlers.items():
                        futures.append(
                            executor.submit(self.process_region, region, vcf_type, vcf)
                        )
                
                # Collect results
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
            # Define headers
            headers = ['chrom', 'start', 'end', 'variant_id', 'position', 'type',
                      'ref', 'alt', 'meth_samples', 'unmeth_samples']
            
            # Write results
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
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--small-vcf', help='Small variants VCF file')
    parser.add_argument('--cnv-vcf', help='Copy number variants VCF file')
    parser.add_argument('--sv-vcf', help='Structural variants VCF file')
    parser.add_argument('--tr-vcf', help='Tandem repeats VCF file')
    args = parser.parse_args()
    
    # Initialize mapper
    mapper = VariantMapper(args.output)
    
    # Add variant files
    vcf_files = {
        'SMALL': args.small_vcf,
        'CNV': args.cnv_vcf,
        'SV': args.sv_vcf,
        'TR': args.tr_vcf
    }
    
    for vcf_type, file_path in vcf_files.items():
        if file_path:
            mapper.add_variant_file(vcf_type, file_path)
    
    if not mapper.vcf_handlers:
        mapper.logger.error("No valid variant files provided")
        exit(1)
    
    # Process variants
    results = mapper.process_bed(args.bed)
    if results:
        if not mapper.write_results(results):
            exit(1)
    else:
        mapper.logger.error("No results generated")
        exit(1)

if __name__ == '__main__':
    main()