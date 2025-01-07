import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def get_region_data(sample_dir):
   data = {}
   for hap in [1, 2]:
       bed_files = list(Path(sample_dir).glob(f"*.hap{hap}.meth_regions.bed"))
       bed_file = bed_files[0]
       df = pd.read_csv(bed_file, sep='\t', comment='#',
                       names=['chrom', 'start', 'end', 'summary_label'])
       df['size'] = df['end'] - df['start']
       
       for label in ['M', 'U']:
           sizes = df[df['summary_label'] == label]['size'].values
           data[f'Hap{hap}_{label}'] = sizes
           data[f'Hap{hap}_{label}_count'] = len(sizes)
   return data

def create_size_plots(data_dir, output_dir):
   samples = sorted([d for d in Path(data_dir).glob('SPM*') if d.is_dir()])
   sample_names = [s.name for s in samples]
   all_data = {sample.name: get_region_data(sample) for sample in samples}

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 6))
   axes = {'M': ax1, 'U': ax2}
   
   colors = {
       'Hap1_M': '#FF0000', 'Hap1_U': '#0000FF',
       'Hap2_M': '#FF6666', 'Hap2_U': '#6666FF'
   }
   
   ylims = {'M': (0, 2050), 'U': (0, 900)}
   x = np.arange(len(sample_names))
   width = 0.35

   for region in ['M', 'U']:
       ax = axes[region]
       
       for i, hap in enumerate(['Hap1', 'Hap2']):
           means = [np.mean(all_data[name][f'{hap}_{region}']) for name in sample_names]
           ax.bar(x + i*width, means, width,
                 color=colors[f'{hap}_{region}'], alpha=0.7,
                 label=f'{hap}')
                 
       ax.set_title(f'Region {region}')
       ax.set_xticks(x + width/2)
       ax.set_xticklabels(sample_names, rotation=45, ha='right')
       ax.set_ylim(ylims[region])
       ax.grid(True, linestyle='--', alpha=0.5, axis='y')
       ax.set_ylabel('Region Size (bp)')
       ax.legend(bbox_to_anchor=(0.5, -0.25), loc='upper center', ncol=2)

   plt.tight_layout()
   plt.savefig(Path(output_dir)/'region_sizes.png', dpi=300, bbox_inches='tight')
   plt.close()

def create_region_counts(data_dir, output_dir):
   samples = sorted([d for d in Path(data_dir).glob('SPM*') if d.is_dir()])
   all_data = {sample.name: get_region_data(sample) for sample in samples}
   df = pd.DataFrame.from_dict(all_data, orient='index')

   fig = plt.figure(figsize=(12, 6))
   x = np.arange(len(df.index))
   width = 0.2
   
   colors = {
       ('Hap1', 'M'): '#FF0000', ('Hap2', 'M'): '#FF6666',
       ('Hap1', 'U'): '#0000FF', ('Hap2', 'U'): '#6666FF'
   }

   for i, (hap, label) in enumerate([
       ('Hap1', 'M'), ('Hap2', 'M'),
       ('Hap1', 'U'), ('Hap2', 'U')
   ]):
       plt.bar(x + i*width, df[f'{hap}_{label}_count'], width,
               label=f'{hap}-{label}', color=colors[(hap, label)])

   plt.xticks(x + 1.5*width, df.index, rotation=45, ha='right')
   plt.legend(bbox_to_anchor=(0.5, -0.25), loc='upper center', ncol=4)
   plt.title('Region Counts by Sample')
   plt.ylabel('Count')
   plt.grid(True, alpha=0.3)
   plt.tight_layout()
   
   plt.savefig(Path(output_dir)/'region_counts.png', dpi=300, bbox_inches='tight')
   plt.close()

def main():
   import sys
   if len(sys.argv) != 3:
       print("Usage: script.py <data_dir> <output_dir>")
       sys.exit(1)
       
   data_dir = sys.argv[1]
   output_dir = sys.argv[2]
   
   Path(output_dir).mkdir(parents=True, exist_ok=True)
   
   create_size_plots(data_dir, output_dir)
   create_region_counts(data_dir, output_dir)

if __name__ == "__main__":
   main()