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

@dataclass
class CpGEffects:
    """Store CpG effects per haplotype"""
    ref_destroyed: Set[int] = field(default_factory=set)  # CpGs destroyed by variants
    denovo_created: Set[int] = field(default_factory=set)  # New CpGs created by variants
    homozygous_denovo: Set[int] = field(default_factory=set)  # New CpGs created by homozygous variants
    phantom_sites: Set[int] = field(default_factory=set)  # Positions that appear as CpGs but aren't in reference
    variant_positions: Set[int] = field(default_factory=set)  # Track actual variant positions

class CpGVariantFilter:
    def __init__(self, vcf_file: str, ref_fasta: str, output_prefix: str, output_dir: str):
        """Initialize filter with input files and output directories"""
        self.vcf = Path(vcf_file)
        self.ref_fasta = Path(ref_fasta)
        self.prefix = output_prefix
        self.setup_directories(output_dir)
        self.setup_logging()
        self.effects = {
            'hap1': CpGEffects(),
            'hap2': CpGEffects()
        }
        self.filter_reasons = defaultdict(lambda: defaultdict(set))

    def setup_directories(self, output_dir: str) -> None:
        """Setup output directory structure"""
        self.output_dir = Path(output_dir)
        self.beds_dir = self.output_dir / 'filtered_beds'
        self.plots_dir = self.output_dir / 'plots'
        self.logs_dir = self.output_dir / 'logs'
        self.passing_dir = self.beds_dir / 'passing'
        self.excluded_dir = self.beds_dir / 'excluded'
        
        for dir in [self.beds_dir, self.plots_dir, self.logs_dir,
                    self.passing_dir, self.excluded_dir]:
            dir.mkdir(parents=True, exist_ok=True)

    def setup_logging(self) -> None:
        """Configure logging with detailed debug level"""
        log_file = self.logs_dir / f'{self.prefix}.filter.detailed.log'
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )

    def is_reference_cpg(self, chrom: str, pos: int, ref_fasta: pysam.FastaFile) -> bool:
        """Check if a position is part of a CpG site in the reference sequence"""
        try:
            # Get 2bp context
            seq = ref_fasta.fetch(chrom, pos, pos + 2).upper()
            return seq == 'CG'
        except Exception as e:
            logging.error(f"Error checking reference CpG at {chrom}:{pos}: {str(e)}")
            return False

    def check_cpg_context(self, chrom: str, pos: int, ref: str, alt: str, 
                         window: int = 2) -> Dict:
        """Analyze CpG context around variant with reference validation"""
        try:
            with pysam.FastaFile(str(self.ref_fasta)) as fasta:
                start = max(0, pos - window)
                end = pos + len(ref) + window
                ref_seq = fasta.fetch(chrom, start, end).upper()
                
                # Create alternate sequence
                var_pos = pos - start
                alt_seq = ref_seq[:var_pos] + alt + ref_seq[var_pos + len(ref):]
                
                # Find CpG positions
                ref_cpgs = self._find_cpg_positions(ref_seq, start)
                alt_cpgs = self._find_cpg_positions(alt_seq, start)
                
                # Validate reference CpGs
                validated_ref_cpgs = {p for p in ref_cpgs 
                                    if self.is_reference_cpg(chrom, p, fasta)}
                
                logging.debug(f"Variant at {chrom}:{pos+1}")
                logging.debug(f"  Reference context: {ref_seq}")
                logging.debug(f"  Alternate context: {alt_seq}")
                logging.debug(f"  Reference CpGs: {sorted(ref_cpgs)}")
                logging.debug(f"  Alternate CpGs: {sorted(alt_cpgs)}")
                logging.debug(f"  Validated ref CpGs: {sorted(validated_ref_cpgs)}")
                
                return {
                    'ref': validated_ref_cpgs,
                    'alt': alt_cpgs,
                    'variant_pos': pos
                }
                
        except Exception as e:
            logging.error(f"Error checking CpG context at {chrom}:{pos}: {str(e)}")
            return None

    @staticmethod
    def _find_cpg_positions(sequence: str, offset: int) -> Set[int]:
        """Find CpG positions in sequence"""
        return {offset + i for i in range(len(sequence)-1) 
                if sequence[i:i+2] == 'CG'}

    def process_variant_effects(self, var, cpg_effects: Dict) -> None:
        """Process variant effects on CpG sites with phantom site detection"""
        destroyed = cpg_effects['ref'] - cpg_effects['alt']
        created = cpg_effects['alt'] - cpg_effects['ref']
        variant_pos = cpg_effects['variant_pos']
        
        if not (destroyed or created):
            return
            
        sample = var.samples[0]
        gt = sample['GT']
        
        logging.debug(f"\nProcessing variant at {var.chrom}:{var.pos}")
        logging.debug(f"  Genotype: {gt}")
        logging.debug(f"  Destroyed CpGs: {sorted(destroyed)}")
        logging.debug(f"  Created CpGs: {sorted(created)}")
        
        if gt == (1,1):  # Homozygous ALT
            if destroyed:
                for hap in ['hap1', 'hap2']:
                    self.effects[hap].ref_destroyed.update(destroyed)
                    self.effects[hap].variant_positions.add(variant_pos)
                    logging.debug(f"  Added destroyed CpGs to {hap}")
            if created:
                for hap in ['hap1', 'hap2']:
                    self.effects[hap].homozygous_denovo.update(created)
                    self.effects[hap].variant_positions.add(variant_pos)
                    logging.debug(f"  Added homozygous denovo CpGs to {hap}")
                    
        elif sample.phased and 'PS' in sample:  # Heterozygous phased
            # Determine which haplotype has the alt and ref alleles
            alt_hap = f'hap{1 if gt[0] == 1 else 2}'
            ref_hap = f'hap{2 if gt[0] == 1 else 1}'
            
            if destroyed:
                self.effects[alt_hap].ref_destroyed.update(destroyed)
                self.effects[alt_hap].variant_positions.add(variant_pos)
                logging.debug(f"  Added destroyed CpGs to {alt_hap}")
                
            if created:
                # For each created CpG site
                for pos in created:
                    # Add as denovo to alt haplotype
                    self.effects[alt_hap].denovo_created.add(pos)
                    self.effects[alt_hap].variant_positions.add(variant_pos)
                    logging.debug(f"  Added created CpG at {pos} to {alt_hap}")
                    
                    # Check if position is actually a CpG in reference
                    with pysam.FastaFile(str(self.ref_fasta)) as fasta:
                        if not self.is_reference_cpg(var.chrom, pos, fasta):
                            # If not a reference CpG, mark as phantom on ref haplotype
                            self.effects[ref_hap].phantom_sites.add(pos)
                            logging.debug(f"  Marked phantom CpG at {pos} on {ref_hap}")
                        else:
                            logging.debug(f"  Position {pos} is real CpG in reference, not marking as phantom")

    def analyze_variants(self) -> None:
        """Analyze variants for CpG effects"""
        vcf = pysam.VariantFile(str(self.vcf))
        var_count = effect_count = 0
        
        for var in vcf.fetch():
            var_count += 1
            if None in var.samples[0]['GT'] or var.samples[0]['GT'] == (0,0):
                continue
                
            pos = var.pos - 1
            alt = var.alts[0] if var.alts else None
            if not alt:
                continue
                
            effects = self.check_cpg_context(var.chrom, pos, var.ref, alt)
            if effects:
                self.process_variant_effects(var, effects)
                effect_count += 1
                
        vcf.close()
        logging.info(f"Processed {var_count:,} variants, found {effect_count:,} affecting CpGs")

    def filter_sites(self, df: pd.DataFrame, haplotype: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Filter methylation sites based on variant effects and phantom sites"""
        exclude_sites = set()
        other_hap = 'hap2' if haplotype == 'hap1' else 'hap1'
        
        # Add sites to exclude
        exclude_sites.update(self.effects[haplotype].ref_destroyed)  # Variants destroying CpGs
        exclude_sites.update(self.effects[haplotype].phantom_sites)  # Phantom CpG sites
        
        # Track reasons
        self.filter_reasons[haplotype]['destroyed'].update(self.effects[haplotype].ref_destroyed)
        self.filter_reasons[haplotype]['phantom'].update(self.effects[haplotype].phantom_sites)
        
        # Filter data
        mask = ~df['start'].isin(exclude_sites)
        return df[mask].copy(), df[~mask].copy()

    def process_bed_file(self, bed_path: Path, bed_type: str) -> Dict:
        """Process individual bed file"""
        df = pd.read_csv(bed_path, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'mod_score', 
                              'haplotype', 'coverage', 'est_mod_count',
                              'est_unmod_count', 'mod_probability'])
        
        if bed_type in ['hap1', 'hap2']:
            kept, excluded = self.filter_sites(df, bed_type)
        else:  # combined
            all_exclusions = set()
            for hap in ['hap1', 'hap2']:
                all_exclusions.update(self.effects[hap].ref_destroyed)
                all_exclusions.update(self.effects[hap].phantom_sites)
            
            mask = ~df['start'].isin(all_exclusions)
            kept, excluded = df[mask], df[~mask]
        
        # Save results
        kept_file = self.passing_dir / f"{self.prefix}.{bed_type}.bed"
        excl_file = self.excluded_dir / f"{self.prefix}.{bed_type}.bed"
        
        kept.to_csv(kept_file, sep='\t', index=False, header=False)
        if len(excluded) > 0:
            excluded.to_csv(excl_file, sep='\t', index=False, header=False)
        
        # Generate plots
        self.generate_plots(excluded, bed_type, excluded=True)
        self.generate_plots(kept, bed_type, excluded=False)
        
        return {
            'total_sites': len(df),
            'kept_sites': len(kept),
            'excluded_sites': len(excluded)
        }

    def generate_plots(self, df: pd.DataFrame, bed_type: str, excluded: bool = False) -> None:
        """Generate distribution plots"""
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
        """Run the complete filtering pipeline"""
        # Find input BED files
        bed_dir = Path(bed_dir)
        if not bed_dir.is_dir():
            raise ValueError(f"BED directory not found: {bed_dir}")
            
        input_beds = {}
        for bed_type in ['combined', 'hap1', 'hap2']:
            bed_path = bed_dir / f"{self.prefix}.{bed_type}.bed"
            if bed_path.exists():
                input_beds[bed_type] = bed_path
                logging.info(f"Found {bed_type} BED file: {bed_path}")
        
        if not input_beds:
            raise ValueError(f"No BED files found with prefix '{self.prefix}'")
        
        # Analyze variants first
        logging.info("Analyzing variants for CpG effects...")
        self.analyze_variants()
        
        # Process each bed file
        results = {}
        for bed_type, bed_path in input_beds.items():
            logging.info(f"\nProcessing {bed_type} data from {bed_path}")
            results[bed_type] = self.process_bed_file(bed_path, bed_type)
            
        return results

    def write_filtering_report(self) -> None:
        """Write detailed filtering report"""
        report_path = self.output_dir / f"{self.prefix}.filtering_report.txt"
        
        with open(report_path, 'w') as f:
            f.write("CpG Methylation Data Filtering Report\n")
            f.write("==================================\n\n")
            
            for haplotype in ['hap1', 'hap2']:
                f.write(f"\n{haplotype.upper()} Details:\n")
                f.write("-" * 20 + "\n")
                f.write(f"Reference CpGs destroyed: {len(self.effects[haplotype].ref_destroyed)}\n")
                f.write(f"Denovo CpGs created: {len(self.effects[haplotype].denovo_created)}\n")
                f.write(f"Homozygous denovo CpGs: {len(self.effects[haplotype].homozygous_denovo)}\n")
                f.write(f"Phantom CpG sites: {len(self.effects[haplotype].phantom_sites)}\n")
                f.write(f"Variant positions: {len(self.effects[haplotype].variant_positions)}\n")
                
                f.write("\nFiltering Decisions:\n")
                reasons = self.filter_reasons[haplotype]
                f.write(f"  Sites excluded due to destruction: {len(reasons['destroyed'])}\n")
                f.write(f"  Sites excluded as phantom CpGs: {len(reasons['phantom'])}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Filter CpG methylation data based on variant effects')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--ref', required=True, help='Reference FASTA file')
    parser.add_argument('--prefix', required=True, help='Sample prefix for BED files')
    parser.add_argument('--bed-dir', required=True, help='Directory containing BED files')
    parser.add_argument('--outdir', default='filtered_results', help='Output directory')
    
    args = parser.parse_args()
    
    try:
        # Initialize and run filter
        filter = CpGVariantFilter(args.vcf, args.ref, args.prefix, args.outdir)
        results = filter.run(args.bed_dir)
        
        # Write filtering report
        filter.write_filtering_report()
        
        print(f"\nFiltering complete. Results in {args.outdir}")
        print(f"See filtering report at {args.outdir}/{args.prefix}.filtering_report.txt")
        
    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()