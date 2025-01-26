# python/segmentation/plot_test.py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def load_segmentation_data(file_path):
    # Load BED file and extract unmethylated region counts
    bed_data = pd.read_csv(file_path, delimiter='\t', header=None)
    unmethylated_counts = [len(region) for index, region in bed_data[3].items()]  # Extract the fourth column as regions
    return unmethylated_counts

def plot_unmethylated_regions(data):
    plt.hist(data, bins=20, edgecolor='black')
    plt.xlabel('Unmethylated Region Counts')
    plt.ylabel('Frequency')
    plt.title('Distribution of Unmethylated Regions')
    plt.show()

def main():
    file_path = 'path/to/your/segmentation/data.bed'  # Replace with actual path to your BED file
    unmethylated_counts = load_segmentation_data(file_path)
    plot_unmethylated_regions(unmethylated_counts)

if __name__ == "__main__":
    main()
