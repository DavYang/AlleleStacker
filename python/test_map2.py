#!/usr/bin/env python3
"""
Parallelized variant mapper for processing all regions in methylation BED files.
Maps variants from VCFs (small, CNV, SV, TR) to methylation regions with full statistics.
"""
import argparse
import subprocess
import pandas as pd
import os
import logging
from typing import Dict, List, Set, Optional
from collections import defaultdict
import gzip
import csv
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
from tqdm import tqdm


class VariantMapper:
    def __init__(self, output_prefix: str, haplotype: str, max_workers: Optional[int] = None):
        self.output_prefix = output_prefix
        self.haplotype = haplotype
        self.vcf_files = {}
        self.vcf_samples = {}
        self.sample_map = {}
        self.stats = defaultdict(lambda: defaultdict(int))
        self.max_workers = max_workers or min(mp.cpu_count(), 20)  # Match SLURM cpus-per-task

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
        self.logger.info(f"Initializing mapper with {self.max_workers} workers")

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

                            # Build sample mapping with SPM276 handling
                            for s in samples:
                                self.sample_map[s.replace('-', '_')] = s
                                self.sample_map[s] = s
                                if s in ('SPM276_1', 'SPM276-1'):
                                    self.sample_map['SPM276_1'] = s
                                    self.sample_map['SPM276-1'] = s
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
                       meth_samples: Set[str], unmeth_samples: Set[str],
                       all_samples: Set[str]) -> Optional[Dict]:
        """Process single variant with haplotype-specific filtering"""
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

            # Check appropriate haplotype allele
            hap_idx = 0 if self.haplotype == 'H1' else 1
            if int(alleles[hap_idx]) == 0:
                continue

            # Categorize sample
            if sample in meth_samples:
                meth_vars.append(f"{sample}:{gt}")
            elif sample in unmeth_samples:
                unmeth_vars.append(f"{sample}:{gt}")
            elif sample in all_samples:
                other_vars.append(f"{sample}:{gt}")

        if not (meth_vars or unmeth_vars or other_vars):
            return None

        # Get variant type
        if vcf_type == 'small':
            if len(fields[3]) == 1 and len(fields[4]) == 1:
                var_type = 'snp'
            elif len(fields[3]) > len(fields[4]):
                var_type = 'deletion'
            else:
                var_type = 'insertion'
        elif vcf_type == 'tr':
            var_type = 'tr'
        elif vcf_type == 'cnv':
            var_type = 'cnv'
        else:
            var_type = 'sv'

        # Update statistics
        self.stats[vcf_type][var_type] += 1
        self.stats[vcf_type]['total'] += 1
        if meth_vars:
            self.stats[vcf_type]['in_meth_samples'] += 1
        if unmeth_vars:
            self.stats[vcf_type]['in_unmeth_samples'] += 1
        if other_vars:
            self.stats[vcf_type]['in_other_samples'] += 1

        return {
            'region_chr': fields[0],
            'region_start': int(fields[1]),
            'region_end': int(fields[1]) + len(fields[3]) - 1,
            'haplotype': self.haplotype,
            'variant_id': fields[2] if fields[2] != '.' else f"{var_type}_{fields[0]}_{fields[1]}",
            'variant_type': var_type,
            'num_samples': len(meth_vars) + len(unmeth_vars) + len(other_vars),
            'num_meth_samples': len(meth_vars),
            'num_unmeth_samples': len(unmeth_vars),
            'num_other_samples': len(other_vars),
            'meth_samples': ','.join(meth_vars) or '.',
            'unmeth_samples': ','.join(unmeth_vars) or '.',
            'other_samples': ','.join(other_vars) or '.'
        }

    def process_region(self, region: pd.Series) -> List[Dict]:
        """Process variants in a region"""
        meth_samples = self._get_samples(region['methylated_samples'])
        unmeth_samples = self._get_samples(region['unmethylated_samples'])
        all_samples = set().union(*self.vcf_samples.values())

        results = []
        for vcf_type, vcf_path in self.vcf_files.items():
            try:
                cmd = ['tabix', vcf_path, f"{region['chrom']}:{region['start']}-{region['end']}"]
                output = subprocess.check_output(cmd, text=True)

                for line in output.splitlines():
                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        continue

                    variant = self.process_variant(fields, vcf_type, meth_samples, unmeth_samples, all_samples)
                    if variant:
                        results.append(variant)

            except subprocess.CalledProcessError:
                pass
            except Exception as e:
                self.logger.error(f"Error processing {vcf_type} variants in {region['chrom']}:{region['start']}: {str(e)}")

        return results

    def write_outputs(self, results: List[Dict]) -> None:
        """Write results to match expected output format"""
        headers = [
            'region_chr', 'region_start', 'region_end', 'haplotype', 'variant_id',
            'variant_type', 'num_samples', 'num_meth_samples', 'num_unmeth_samples',
            'num_other_samples', 'meth_samples', 'unmeth_samples', 'other_samples'
        ]

        # Write main results file
        output_file = f"{self.output_prefix}_variants.tsv"
        with open(output_file, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
            writer.writeheader()
            for r in results:
                writer.writerow(r)

        # Write statistics
        stats_file = f"{self.output_prefix}_stats.txt"
        with open(stats_file, 'w') as f:
            f.write("Variant Mapping Statistics\n")
            f.write("=========================\n\n")
            
            for vcf_type, type_stats in self.stats.items():
                f.write(f"\n{vcf_type.upper()} Variants:\n")
                f.write("-" * (len(vcf_type) + 10) + "\n")
                f.write(f"Total variants: {type_stats['total']}\n\n")
                
                f.write("By type:\n")
                for var_type in ['snp', 'insertion', 'deletion', 'cnv', 'sv', 'tr']:
                    if type_stats[var_type]:
                        f.write(f"  {var_type}: {type_stats[var_type]}\n")
                
                f.write("\nBy sample category:\n")
                f.write(f"  In methylated samples: {type_stats['in_meth_samples']}\n")
                f.write(f"  In unmethylated samples: {type_stats['in_unmeth_samples']}\n")
                f.write(f"  In other samples: {type_stats['in_other_samples']}\n\n")

    def run(self, bed_path: str) -> bool:
        """Process all regions with progress tracking"""
        try:
            # Read BED file
            df = pd.read_csv(bed_path, sep='\t')
            self.logger.info(f"Processing {len(df)} regions for haplotype {self.haplotype}")

            # Process regions in parallel with progress bar
            results = []
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = []
                
                # Submit all regions for processing
                for _, region in df.iterrows():
                    futures.append(executor.submit(self.process_region, region))
                
                # Process results as they complete
                for future in tqdm(futures, desc=f"Processing {self.haplotype} regions"):
                    try:
                        region_results = future.result()
                        results.extend(region_results)
                    except Exception as e:
                        self.logger.error(f"Error processing region: {e}")

            if not results:
                self.logger.error("No variants found")
                return False

            # Write outputs
            self.write_outputs(results)
            self.logger.info(f"Processed {len(results)} variants")
            return True

        except Exception as e:
            self.logger.error(f"Error running variant mapping: {e}")
            return False


def main():
    parser = argparse.ArgumentParser(description='Map variants to methylation regions')
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
    parser.add_argument('--tr-vcf',
                       help='Tandem repeats VCF')
    parser.add_argument('--threads', type=int,
                       help='Number of processing threads')

    args = parser.parse_args()

    # Run variant mapping
    mapper = VariantMapper(args.output_prefix, args.haplotype, args.threads)
    if not mapper.load_vcfs({
        'small': args.small_vcf,
        'cnv': args.cnv_vcf,
        'sv': args.sv_vcf,
        'tr': args.tr_vcf
    }):
        mapper.logger.error("No valid VCF files provided")
        exit(1)

    if not mapper.run(args.bed):
        mapper.logger.error("Variant mapping failed")
        exit(1)


if __name__ == '__main__':
    main()