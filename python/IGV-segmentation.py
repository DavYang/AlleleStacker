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
        'other': '0,255,0'  # Green for others
    }

    # Process each meth_regions.bed file
    for filename in os.listdir(input_dir):
        if filename.endswith('meth_regions.bed'):
            file_path = os.path.join(input_dir, filename)
            df = pd.read_csv(file_path, sep='\t', comment='#', header=None, 
                           names=['chrom', 'start', 'end', 'summary_label'])

            # Split into H1 and H2 dataframes
            h1_df = df[df['summary_label'].str.startswith('H1_')].copy()
            h2_df = df[df['summary_label'].str.startswith('H2_')].copy()

            # Function to sort and format dataframe
            def prepare_output_df(df):
                # Sort chromosomes naturally
                df['chrom_num'] = df['chrom'].str.extract('(\d+|X|Y|M)').fillna('0')
                df['chrom_num'] = pd.Categorical(df['chrom_num'], 
                                               categories=[str(i) for i in range(23)] + ['X', 'Y', 'M'], 
                                               ordered=True)
                
                # Sort by chromosome and position
                df = df.sort_values(['chrom_num', 'start', 'end'])
                
                # Create output dataframe with exact structure
                return pd.DataFrame({
                    'chrom': df['chrom'],
                    'start': df['start'],
                    'end': df['end'],
                    'name': df['summary_label'],
                    'score': 0,
                    'strand': '.',
                    'thickStart': df['start'],
                    'thickEnd': df['end'],
                    'rgb': df['summary_label'].apply(
                        lambda x: COLOR_MAP['M'] if x.endswith('M') 
                                 else (COLOR_MAP['U'] if x.endswith('U') 
                                 else COLOR_MAP['other'])
                    )
                })

            # Process and write H1 regions
            if not h1_df.empty:
                h1_output = prepare_output_df(h1_df)
                h1_output_path = os.path.join(output_dir, f'H1_{filename}')
                h1_output.to_csv(h1_output_path, sep='\t', index=False, header=False)

            # Process and write H2 regions
            if not h2_df.empty:
                h2_output = prepare_output_df(h2_df)
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