#!/usr/bin/env python
import argparse
import pysam 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import logging
import sys
from dataclasses import dataclass, field
from typing import Set, Dict, List, Tuple
from collections import defaultdict
import multiprocessing as mp
from functools import partial

def create_defaultdict():
    """Helper function to create defaultdict for filter reasons"""
    return defaultdict(set)

@dataclass
class CpGEffects:
    """Store CpG effects per haplotype"""
    ref_destroyed: Set[int] = field(default_factory=set)
    denovo_created: Set[int] = field(default_factory=set)
    homozygous_denovo: Set[int] = field(default_factory=set)
    phantom_sites: Set[int] = field(default_factory=set)
    variant_positions: Set[int] = field(default_factory=set)

@dataclass
class VariantData:
    """Serializable variant data class"""
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: tuple
    phased: bool
    ps: bool

def convert_variant(var: pysam.VariantRecord) -> VariantData:
    """Convert pysam variant to serializable format"""
    return VariantData(
        chrom=var.chrom,
        pos=var.pos,
        ref=var.ref,
        alt=var.alts[0] if var.alts else None,
        gt=var.samples[0]['GT'],
        phased=var.samples[0].phased,
        ps='PS' in var.samples[0]
    )

def process_variant_chunk(chunk_data):
    """Process a chunk of variants - designed for multiprocessing"""
    ref_fasta_path, variant_data_list = chunk_data
    
    local_effects = {
        'hap1': CpGEffects(),
        'hap2': CpGEffects()
    }
    
    with pysam.FastaFile(str(ref_fasta_path)) as fasta:
        for var in variant_data_list:
            if None in var.gt or var.gt == (0,0):
                continue
                
            if not var.alt:
                continue
            
            try:
                # Check CpG context
                start = max(0, var.pos - 1 - 2)
                end = var.pos - 1 + len(var.ref) + 2
                ref_seq = fasta.fetch(var.chrom, start, end).upper()
                
                var_pos = 2  # Window size is 2
                alt_seq = ref_seq[:var_pos] + var.alt + ref_seq[var_pos + len(var.ref):]
                
                # Find CpG positions
                ref_cpgs = set()
                alt_cpgs = set()
                
                for i in range(len(ref_seq)-1):
                    if ref_seq[i:i+2] == 'CG':
                        ref_cpgs.add(start + i)
                    if alt_seq[i:i+2] == 'CG':
                        alt_cpgs.add(start + i)
                
                # Validate reference CpGs
                validated_ref_cpgs = set()
                for p in ref_cpgs:
                    try:
                        if fasta.fetch(var.chrom, p, p+2).upper() == 'CG':
                            validated_ref_cpgs.add(p)
                    except:
                        continue
                
                # Process effects
                destroyed = validated_ref_cpgs - alt_cpgs
                created = alt_cpgs - validated_ref_cpgs
                variant_pos = var.pos - 1
                
                if var.gt == (1,1):  # Homozygous ALT
                    if destroyed:
                        for hap in ['hap1', 'hap2']:
                            local_effects[hap].ref_destroyed.update(destroyed)
                            local_effects[hap].variant_positions.add(variant_pos)
                    if created:
                        for hap in ['hap1', 'hap2']:
                            local_effects[hap].homozygous_denovo.update(created)
                            local_effects[hap].variant_positions.add(variant_pos)
                            
                elif var.phased and var.ps:  # Heterozygous phased
                    alt_hap = f'hap{1 if var.gt[0] == 1 else 2}'
                    ref_hap = f'hap{2 if var.gt[0] == 1 else 1}'
                    
                    if destroyed:
                        local_effects[alt_hap].ref_destroyed.update(destroyed)
                        local_effects[alt_hap].variant_positions.add(variant_pos)
                        
                    if created:
                        for pos in created:
                            local_effects[alt_hap].denovo_created.add(pos)
                            local_effects[alt_hap].variant_positions.add(variant_pos)
                            
                            try:
                                if fasta.fetch(var.chrom, pos, pos+2).upper() != 'CG':
                                    local_effects[ref_hap].phantom_sites.add(pos)
                            except:
                                continue
                
            except Exception as e:
                logging.error(f"Error processing variant at {var.chrom}:{var.pos}: {str(e)}")
                continue
    
    return local_effects

class CpGVariantFilter:
    def __init__(self, vcf_file: str, ref_fasta: str, output_prefix: str, output_dir: str):
        self.vcf = Path(vcf_file)
        self.ref_fasta = Path(ref_fasta)
        self.prefix = output_prefix
        self.setup_directories(output_dir)
        self._setup_logging()
        self.effects = {
            'hap1': CpGEffects(),
            'hap2': CpGEffects()
        }
        self.results = defaultdict(dict)  # Store results for each bed type
        self.filter_reasons = defaultdict(create_defaultdict)
        self.preserved_sites = defaultdict(set)  # Track preserved high-confidence sites

    def setup_directories(self, output_dir: str) -> None:
        self.output_dir = Path(output_dir)
        self.beds_dir = self.output_dir / 'filtered_beds'
        self.plots_dir = self.output_dir / 'plots'
        self.logs_dir = self.output_dir / 'logs'
        self.passing_dir = self.beds_dir / 'passing'
        self.excluded_dir = self.beds_dir / 'excluded'
        
        for dir in [self.beds_dir, self.plots_dir, self.logs_dir,
                    self.passing_dir, self.excluded_dir]:
            dir.mkdir(parents=True, exist_ok=True)

    def _setup_logging(self) -> None:
        log_file = self.logs_dir / f'{self.prefix}.filter.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        logging.info(f"Starting CpG filtering for {self.prefix}")
        logging.info(f"VCF file: {self.vcf}")
        logging.info(f"Reference: {self.ref_fasta}")

    def analyze_variants(self) -> None:
        variants = []
        vcf = pysam.VariantFile(str(self.vcf))
        for var in vcf.fetch():
            variants.append(convert_variant(var))
        vcf.close()
        
        chunk_size = max(1000, len(variants) // (mp.cpu_count() * 2))
        variant_chunks = [variants[i:i + chunk_size] for i in range(0, len(variants), chunk_size)]
        chunk_data = [(str(self.ref_fasta), chunk) for chunk in variant_chunks]
        
        with mp.Pool() as pool:
            chunk_results = pool.map(process_variant_chunk, chunk_data)
        
        for chunk_effect in chunk_results:
            for hap in ['hap1', 'hap2']:
                for attr in ['ref_destroyed', 'denovo_created', 'homozygous_denovo', 
                           'phantom_sites', 'variant_positions']:
                    getattr(self.effects[hap], attr).update(
                        getattr(chunk_effect[hap], attr))

    def filter_sites(self, df: pd.DataFrame, haplotype: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Filter methylation sites based on variant effects and phantom sites"""
        exclude_sites = set()
        
        # Add sites to exclude
        exclude_sites.update(self.effects[haplotype].ref_destroyed)
        exclude_sites.update(self.effects[haplotype].phantom_sites)
        
        # Track reasons
        self.filter_reasons[haplotype]['destroyed'].update(self.effects[haplotype].ref_destroyed)
        self.filter_reasons[haplotype]['phantom'].update(self.effects[haplotype].phantom_sites)

        # Track preserved high-confidence sites
        high_conf_mask = (df['mod_probability'] > 0.9) & (df['coverage'] >= 10)
        preserve_mask = high_conf_mask & ~df['start'].isin(self.effects[haplotype].variant_positions)
        preserved_sites = set(df[preserve_mask & df['start'].isin(exclude_sites)]['start'])
        self.preserved_sites[haplotype].update(preserved_sites)
        
        # Final filtering
        exclude_mask = df['start'].isin(exclude_sites)
        final_mask = ~exclude_mask | preserve_mask
        
        return df[final_mask].copy(), df[~final_mask].copy()

    def process_bed_file(self, bed_path: Path, bed_type: str) -> Dict:
        chunk_size = 100000
        chunks = pd.read_csv(bed_path, sep='\t', header=None, chunksize=chunk_size,
                           names=['chrom', 'start', 'end', 'mod_score', 
                                 'haplotype', 'coverage', 'est_mod_count',
                                 'est_unmod_count', 'mod_probability'])
        
        kept_chunks = []
        excluded_chunks = []
        
        for chunk in chunks:
            if bed_type in ['hap1', 'hap2']:
                kept, excluded = self.filter_sites(chunk, bed_type)
            else:
                all_exclusions = set()
                for hap in ['hap1', 'hap2']:
                    all_exclusions.update(self.effects[hap].ref_destroyed)
                    all_exclusions.update(self.effects[hap].phantom_sites)
                
                mask = ~chunk['start'].isin(all_exclusions)
                kept, excluded = chunk[mask], chunk[~mask]
            
            kept_chunks.append(kept)
            excluded_chunks.append(excluded)
        
        kept = pd.concat(kept_chunks, ignore_index=True)
        excluded = pd.concat(excluded_chunks, ignore_index=True)
        
        stats = {
            'total_sites': len(kept) + len(excluded),
            'kept_sites': len(kept),
            'excluded_sites': len(excluded),
        }
        
        if bed_type in ['hap1', 'hap2']:
            stats.update({
                'destroyed_sites': len(self.effects[bed_type].ref_destroyed),
                'phantom_sites': len(self.effects[bed_type].phantom_sites),
                'preserved_sites': len(self.preserved_sites[bed_type])
            })
        
        # Store results
        self.results[bed_type] = stats
        
        # Save filtered files
        kept_file = self.passing_dir / f"{self.prefix}.{bed_type}.filtered.bed"
        excl_file = self.excluded_dir / f"{self.prefix}.{bed_type}.excluded.bed"
        
        kept.to_csv(kept_file, sep='\t', index=False, header=False)
        if not excluded.empty:
            excluded.to_csv(excl_file, sep='\t', index=False, header=False)
            self.generate_plots(excluded, bed_type, excluded=True)
        
        self.generate_plots(kept, bed_type, excluded=False)
        
        return stats

    def write_summary_csv(self) -> None:
        """Write summary statistics to CSV file"""
        summary_file = self.output_dir / "filtering_summary.csv"
        
        # Create summary data for this sample
        summary_data = {
            'Sample': self.prefix,
            'Total CpG Sites (Combined)': self.results['combined']['total_sites'],
            'Kept CpG Sites (Combined)': self.results['combined']['kept_sites'],
            'Excluded CpG Sites (Combined)': self.results['combined']['excluded_sites']
        }
        
        # Add haplotype-specific data
        for hap in ['hap1', 'hap2']:
            if hap in self.results:
                summary_data.update({
                    f'Total CpG Sites ({hap})': self.results[hap]['total_sites'],
                    f'Kept CpG Sites ({hap})': self.results[hap]['kept_sites'],
                    f'Excluded CpG Sites ({hap})': self.results[hap]['excluded_sites'],
                    f'Destroyed CpGs ({hap})': len(self.effects[hap].ref_destroyed),
                    f'Phantom CpGs ({hap})': len(self.effects[hap].phantom_sites),
                    f'High Confidence Preserved ({hap})': len(self.preserved_sites[hap])
                })
        
        # Create DataFrame and save
        summary_df = pd.DataFrame([summary_data])
        
        # If file exists, append; otherwise create new
        if summary_file.exists():
            try:
                existing_df = pd.read_csv(summary_file)
                summary_df = pd.concat([existing_df, summary_df], ignore_index=True)
            except pd.errors.EmptyDataError:
                pass
        
        # Sort by sample name and save
        summary_df = summary_df.sort_values('Sample')
        summary_df.to_csv(summary_file, index=False)
        logging.info(f"Written summary statistics to {summary_file}")

    def generate_plots(self, df: pd.DataFrame, bed_type: str, excluded: bool = False) -> None:
        if len(df) == 0:
            return
            
        plot_type = 'excluded' if excluded else 'passing'
        for metric in ['mod_probability', 'coverage']:
            plt.figure(figsize=(10, 6))
            sns.histplot(data=df[metric], bins=50)
            plt.title(f'{metric.replace("_", " ").title()} Distribution - {bed_type} ({plot_type})')
            plt.xlabel(metric.replace("_", " ").title())
            plt.ylabel('Count')
            plt.savefig(self.plots_dir / f'{self.prefix}.{bed_type}.{plot_type}.{metric}.pdf')
            plt.close()

    def run(self, bed_dir: str) -> Dict:
        """Run optimized pipeline"""
        bed_dir = Path(bed_dir)
        if not bed_dir.is_dir():
            raise ValueError(f"BED directory not found: {bed_dir}")
            
        input_beds = {}
        for bed_type in ['combined', 'hap1', 'hap2']:
            bed_path = bed_dir / f"{self.prefix}.{bed_type}.bed"
            if bed_path.exists():
                input_beds[bed_type] = bed_path
                logging.info(f"Found {bed_type} BED file: {bed_path}")
            else:
                logging.warning(f"BED file not found: {bed_path}")
        
        if not input_beds:
            raise ValueError(f"No BED files found with prefix '{self.prefix}'")
        
        logging.info("Analyzing variants...")
        self.analyze_variants()
        
        logging.info("Processing bed files...")
        results = {}
        for bed_type, bed_path in input_beds.items():
            logging.info(f"Processing {bed_type} data...")
            results[bed_type] = self.process_bed_file(bed_path, bed_type)
        
        # Write summary CSV with all results
        self.write_summary_csv()
        
        return results

def main():
    parser = argparse.ArgumentParser(
        description='Filter CpG methylation data based on variant effects')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--ref', required=True, help='Reference FASTA file')
    parser.add_argument('--prefix', required=True, help='Sample prefix for BED files')
    parser.add_argument('--bed-dir', required=True, help='Directory containing BED files')
    parser.add_argument('--outdir', default='filtered_results', help='Output directory')
    parser.add_argument('--threads', type=int, default=mp.cpu_count(),
                       help='Number of threads for parallel processing')
    
    args = parser.parse_args()
    
    try:
        # Set number of available CPU cores for parallel processing
        if args.threads > 0:
            mp.set_start_method('spawn', force=True)
        
        # Initialize and run filter
        filter = CpGVariantFilter(args.vcf, args.ref, args.prefix, args.outdir)
        results = filter.run(args.bed_dir)
        filter.write_filtering_report()
        
        logging.info("Processing complete")
        print(f"\nFiltering complete. Results in {args.outdir}")
        print(f"See filtering report at {args.outdir}/{args.prefix}.filtering_report.txt")
        print(f"See summary CSV at {args.outdir}/filtering_summary.csv")
        
    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()