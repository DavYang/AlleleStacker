# segmentation_group_summary.py
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)

def generate_summary(sample_dir, output_file):
   samples = sorted([d for d in Path(sample_dir).glob('SPM*') if d.is_dir()])
   data = []
   
   for sample_path in samples:
       sample_name = sample_path.name
       hap1_file = sample_path / f'{sample_name}.hap1.meth_regions.bed'
       hap2_file = sample_path / f'{sample_name}.hap2.meth_regions.bed'
       
       if hap1_file.exists() and hap2_file.exists():
           logging.info(f"Processing {sample_name}")
           
           row = {'sample': sample_name}
           
           for hap, file in [('Hap1', hap1_file), ('Hap2', hap2_file)]:
               regions = pd.read_csv(file, sep='\t')
               for label in ['M', 'U']:
                   subset = regions[regions['summary_label'] == label]
                   size = subset['end'] - subset['start']
                   
                   prefix = f"{hap}_{label}"
                   row[f'{prefix}_count'] = len(subset)
                   row[f'{prefix}_mean_size'] = size.mean() if len(size) > 0 else 0
                   row[f'{prefix}_std_size'] = size.std() if len(size) > 0 else 0
           
           data.append(row)
   
   df = pd.DataFrame(data)
   df.set_index('sample').to_csv(output_file)
   logging.info(f"Summary saved to {output_file}")

if __name__ == "__main__":
   import sys
   if len(sys.argv) != 3:
       sys.exit("Usage: python segmentation_group_summary.py <sample_dir> <output_file>")
   
   generate_summary(sys.argv[1], sys.argv[2])