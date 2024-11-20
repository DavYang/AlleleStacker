import os
import pandas as pd
import sys

def filter_meth_regions(input_dir, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define IGV-compatible RGB colors
    COLOR_MAP = {
        'M': '255,0,0',    # Red for methylated
        'U': '0,0,255',    # Blue for unmethylated
        'other': '0,255,0'  # Green for others, matching your example
    }

    # Process each meth_regions.bed file
    for filename in os.listdir(input_dir):
        if filename.endswith('meth_regions.bed'):
            file_path = os.path.join(input_dir, filename)
            df = pd.read_csv(file_path, sep='\t', comment='#', header=None, 
                           names=['chrom', 'start', 'end', 'summary_label'])
            
            # Filter and process H1 regions
            h1_df = df[df['summary_label'].str.startswith('H1_')]
            if not h1_df.empty:
                h1_output = pd.DataFrame({
                    'chrom': h1_df['chrom'],
                    'start': h1_df['start'],
                    'end': h1_df['end'],
                    'name': h1_df['summary_label'],
                    'score': 0,  # Set score to 0 as in example
                    'strand': '.',
                    'thickStart': h1_df['start'],
                    'thickEnd': h1_df['end'],
                    'rgb': h1_df['summary_label'].apply(
                        lambda x: COLOR_MAP['M'] if x.endswith('M') 
                                 else (COLOR_MAP['U'] if x.endswith('U') 
                                 else COLOR_MAP['other'])
                    )
                })
                h1_output_path = os.path.join(output_dir, f'H1_{filename}')
                h1_output.to_csv(h1_output_path, sep='\t', index=False, header=False)

            # Filter and process H2 regions
            h2_df = df[df['summary_label'].str.startswith('H2_')]
            if not h2_df.empty:
                h2_output = pd.DataFrame({
                    'chrom': h2_df['chrom'],
                    'start': h2_df['start'],
                    'end': h2_df['end'],
                    'name': h2_df['summary_label'],
                    'score': 0,  # Set score to 0 as in example
                    'strand': '.',
                    'thickStart': h2_df['start'],
                    'thickEnd': h2_df['end'],
                    'rgb': h2_df['summary_label'].apply(
                        lambda x: COLOR_MAP['M'] if x.endswith('M') 
                                 else (COLOR_MAP['U'] if x.endswith('U') 
                                 else COLOR_MAP['other'])
                    )
                })
                h2_output_path = os.path.join(output_dir, f'H2_{filename}')
                h2_output.to_csv(h2_output_path, sep='\t', index=False, header=False)

    print("Filtering complete. Check the output directory for results.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_meth_regions.py <input_dir> <output_dir>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    filter_meth_regions(input_directory, output_directory)