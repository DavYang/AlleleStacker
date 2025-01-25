import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import logging

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.getLogger().addHandler(logging.StreamHandler())

def load_methylation_data(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
        required_cols = ['chrom', 'start', 'end']
        if not all(col in df.columns for col in required_cols):
            logging.error(f"Missing required columns in {file_path}")
            return None
        df['size'] = df['end'] - df['start']
        return df
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return None

def check_input_dir(input_dir):
    if not Path(input_dir).is_dir():
        raise FileNotFoundError(f"Input directory '{input_dir}' does not exist.")

def check_output_dir(output_dir):
    if not Path(output_dir).is_dir():
        raise FileNotFoundError(f"Output directory '{output_dir}' does not exist.")

def main():
    parser = argparse.ArgumentParser(description="Plot methylation regions by size and position.")
    parser.add_argument("--input-dir", required=True, help="Path to input directory")
    parser.add_argument("--output-dir", required=True, help="Path to output directory")
    parser.add_argument("--sample-name", required=True, help="Sample name to process")

    args = parser.parse_args()

    check_input_dir(args.input_dir)
    check_output_dir(args.output_dir)

    input_files = list(Path(args.input_dir).glob("*.tsv"))
    if not input_files:
        raise FileNotFoundError(f"No TSV files found in '{args.input_dir}'.")

    for file in input_files:
        df = load_methylation_data(file)
        if df is None:
            continue

        sample_name = Path(file).stem
        if sample_name != args.sample_name:
            logging.warning(f"Skipping file '{file}' as it does not belong to sample '{args.sample_name}'.")
            continue

        plot_distribution(df, args.output_dir)

def plot_distribution(df, output_dir):
    try:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        figure = plt.figure()
        sns.boxplot(x="size", y="start", data=df)
        figure.savefig(Path(output_dir)/f"{args.sample_name}_size_position.png")
        logging.info(f"Saved size-position plot for sample '{args.sample_name}' to {Path(output_dir)/f'{args.sample_name}_size_position.png'}.")
    except Exception as e:
        logging.error(f"Error creating or saving figure: {e}")

if __name__ == "__main__":
    main()
