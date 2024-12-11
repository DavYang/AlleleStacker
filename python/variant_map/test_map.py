#!/usr/bin/env python3
"""
test_variant_mapper.py - Maps variants to methylation regions with simplified output and haplotype-specific filtering.
"""
import argparse
import subprocess
import pandas as pd
import os
import logging
from typing import Dict, List, Set, Optional
from collections import defaultdict
import gzip
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
import csv


class TestMapper:
    def __init__(self, output_prefix: str, haplotype: str):
        self.output_prefix = output_prefix
        self.haplotype = haplotype
        self.vcf_files = {}
        self.vcf_samples = {}
        self.sample_map = {}
        self.stats = {
            'snp': 0,
            'insertion': 0,
            'deletion': 0,
            'cnv': 0,
            'sv': 0
        }

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

    def load_vcfs(self, vcf_files: Dict[str, str]) -> bool:
        """Load and validate VCF files"""
        for vcf_type, path in vcf_files.items():
            if not path or not all(os.path.exists(f) for f in [path, f"{path}.tbi"]):
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
                            break
                self.logger.info(f"Loaded {vcf_type} VCF with {len(samples)} samples")
            except Exception as e:
                self.logger.error(f"Error loading {vcf_type} VCF: {e}")

        return bool(self.vcf_files)

    def _get_samples(self, sample_str: str) -> Set[str]:
        """Parse and normalize sample names"""
        if pd.isna(sample_str) or sample_str == '.':
            return set()
        return {self.sample_map.get(s.strip(), s.strip())
                for s in sample_str.replace(';', ',').split(',') if s.strip()}

    def process_variant(self, fields: List[str], vcf_type: str,
                        meth_samples: Set[str], unmeth_samples: Set[str]) -> Optional[Dict]:
        """Process single variant"""
        genotypes = dict(zip(self.vcf_samples[vcf_type], fields[9:]))
        meth_vars, unmeth_vars, other_vars = [], [], []

        for sample, gt_str in genotypes.items():
            if gt_str in {'./.', '.|.'}:
                continue
            gt = gt_str.split(':')[0]
            if '|' not in gt:
                continue
            alleles = gt.split('|')
            if len(alleles) != 2:
                continue

            if self.haplotype == 'H1' and int(alleles[0]) != 0:
                if sample in meth_samples:
                    meth_vars.append(f"{sample}:{gt}")
                elif sample in unmeth_samples:
                    unmeth_vars.append(f"{sample}:{gt}")
                else:
                    other_vars.append(f"{sample}:{gt}")
            elif self.haplotype == 'H2' and int(alleles[1]) != 0:
                if sample in meth_samples:
                    meth_vars.append(f"{sample}:{gt}")
                elif sample in unmeth_samples:
                    unmeth_vars.append(f"{sample}:{gt}")
                else:
                    other_vars.append(f"{sample}:{gt}")

        if not (meth_vars or unmeth_vars or other_vars):
            self.logger.debug(f"No haplotype-specific variants found in region {fields[0]}:{fields[1]}")
            return None

        # Get variant type and update stats
        if vcf_type == 'small':
            if len(fields[3]) == 1 and len(fields[4]) == 1:
                var_type = 'snp'
            elif len(fields[3]) > len(fields[4]):
                var_type = 'deletion'
            else:
                var_type = 'insertion'
        elif vcf_type == 'cnv':
            var_type = 'cnv'
        else:
            var_type = 'sv'
        self.stats[var_type] += 1

        return {
            'region_chr': fields[0],
            'region_start': int(fields[1]),
            'region_end': int(fields[1]) + len(fields[3]) - 1,
            'haplotype': self.haplotype,
            'variant_id': fields[2] if fields[2] != '.' else f"{var_type}_{fields[0]}_{fields[1]}",
            'variant_type': var_type,
            'meth_samples': ','.join(meth_vars) or '0',
            'unmeth_samples': ','.join(unmeth_vars) or '0'
        }

    def process_region(self, region: pd.Series) -> List[Dict]:
        # Get methylated and unmethylated samples
        meth_samples = self._get_samples(region['methylated_samples'])
        unmeth_samples = self._get_samples(region['unmethylated_samples'])

        # Debugging: Log sample information from the BED file
        self.logger.debug(f"Region: {region['chrom']}:{region['start']}-{region['end']}")
        self.logger.debug(f"Methylated samples in BED file: {meth_samples}")
        self.logger.debug(f"Unmethylated samples in BED file: {unmeth_samples}")

        results = []

        # Process variants for the region from each VCF type
        for vcf_type, vcf_path in self.vcf_files.items():
            try:
                # Extract variants from VCF using tabix
                cmd = ['tabix', vcf_path, f"{region['chrom']}:{region['start']}-{region['end']}"]
                output = subprocess.check_output(cmd, text=True)

                for line in output.splitlines():
                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        continue

                    # Process the variant
                    variant = self.process_variant(fields, vcf_type, meth_samples, unmeth_samples)
                    if variant:
                        # Debugging: Log sample information for the variant
                        self.logger.debug(
                            f"Variant {variant['variant_id']} in region "
                            f"{region['chrom']}:{region['start']}-{region['end']}: "
                            f"Methylated samples mapped: {variant['meth_samples']}, "
                            f"Unmethylated samples mapped: {variant['unmeth_samples']}"
                        )

                        # Update the variant with region-specific information
                        variant.update({
                            'region_chr': region['chrom'],
                            'region_start': region['start'],
                            'region_end': region['end']
                        })
                        results.append(variant)

            except subprocess.CalledProcessError:
                # Log if no variants found for this VCF type
                self.logger.debug(
                    f"No variants found in {vcf_type} VCF for region "
                    f"{region['chrom']}:{region['start']}-{region['end']}"
                )

        # If this is a single region run (e.g., during testing), log results summary
        if len(results) > 0:
            self.logger.info(f"Processed {len(results)} variants for region "
                             f"{region['chrom']}:{region['start']}-{region['end']}")

        return results


    def write_outputs(self, results: List[Dict]) -> None:
        """Write combined output file"""
        headers = [
            'region_chr', 'region_start', 'region_end',
            'haplotype', 'variant_id', 'variant_type', 'meth_samples', 'unmeth_samples'
        ]

        with open(f"{self.output_prefix}", 'w') as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
            writer.writeheader()
            for r in results:
                writer.writerow(r)

    def run_test(self, bed_path: str, sample_size: int = 50) -> bool:
        """Run the test"""
        try:
            # Read and sample regions
            df = pd.read_csv(bed_path, sep='\t')
            test_regions = df.sample(n=min(sample_size, len(df)))
            self.logger.info(f"Testing {len(test_regions)} regions for haplotype {self.haplotype}")

            # Process regions in parallel
            results = []
            with ThreadPoolExecutor(max_workers=min(mp.cpu_count(), 8)) as executor:
                futures = [executor.submit(self.process_region, region)
                           for _, region in test_regions.iterrows()]

                for f in futures:
                    results.extend(f.result())

            if not results:
                self.logger.error("No variants found")
                return False

            # Write results
            self.write_outputs(results)
            self.logger.info(f"Processed {len(results)} variants")
            return True

        except Exception as e:
            self.logger.error(f"Error running test: {e}")
            return False


def main():
    parser = argparse.ArgumentParser(
        description='Test variant mapping against methylation regions'
    )
    parser.add_argument('--bed', required=True,
                        help='Input BED file with methylation regions')
    parser.add_argument('--haplotype', required=True, choices=['H1', 'H2'],
                        help='Haplotype to process')
    parser.add_argument('--output-prefix', required=True,
                        help='Prefix for output files')
    parser.add_argument('--small-vcf',
                        help='Small variants VCF')
    parser.add_argument('--cnv-vcf',
                        help='CNV variants VCF')
    parser.add_argument('--sv-vcf',
                        help='Structural variants VCF')
    parser.add_argument('--sample-size', type=int, default=50,
                        help='Number of regions to test (default: 50)')

    args = parser.parse_args()

    # Run test
    mapper = TestMapper(args.output_prefix, args.haplotype)
    if not mapper.load_vcfs({
        'small': args.small_vcf,
        'cnv': args.cnv_vcf,
        'sv': args.sv_vcf
    }):
        mapper.logger.error("No valid VCF files provided")
        exit(1)

    if not mapper.run_test(args.bed, args.sample_size):
        mapper.logger.error("Test failed")
        exit(1)


if __name__ == '__main__':
    main()
