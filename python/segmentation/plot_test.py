import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import logging

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.getLogger().addHandler(logging.StreamHandler())

def load_methylation_data(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
        required_cols = ['chrom', 'start', 'end']
        if not all(col in df.columns for col in required_cols):
            logging.error(f"Missing required columns in {file_path}")
            return None
        df['size'] = df['end'] - df['start']
        logging.debug(f"Successfully loaded and processed data from {file_path}")
        return df
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return None

def check_input_dir(input_dir):
    if not Path(input_dir).is_dir():
        raise FileNotFoundError(f"Input directory '{input_dir}' does not exist.")
    logging.debug(f"Confirmed that input directory '{input_dir}' exists")

def check_output_dir(output_dir):
    if not Path(output_dir).is_dir():
        raise FileNotFoundError(f"Output directory '{output_dir}' does not exist.")
    logging.debug(f"Confirmed that output directory '{output_dir}' exists")

def main():
    # rest of the code...

if __name__ == "__main__":
    setup_logging()
    main()
