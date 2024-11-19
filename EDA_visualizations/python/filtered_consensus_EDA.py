import matplotlib.pyplot as plt
from collections import Counter

def parse_bed_file(file_path):
    sample_counts = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                sample_count = int(fields[3])  # Directly use the fourth column as the sample count
                sample_counts.append(sample_count)
    return sample_counts

def plot_frequency_distribution(sample_counts):
    counter = Counter(sample_counts)
    max_count = max(counter.keys())
    
    x = list(range(1, max_count + 1))
    y = [counter[i] for i in x]
    
    #plt.figure(figsize=(15, 8))
    plt.bar(x, y, color='red')  
    plt.xlabel('Number of Samples per Region')
    plt.ylabel('Count')
    # plt.title('Frequency Distribution of Samples per Region (Haplotype 1)')
    plt.title('Frequency Distribution of Samples per Region (Haplotype 2)')
    plt.xticks(range(0, max_count + 1, 1))  # Show all x-axis labels
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    for i, v in enumerate(y):
        if v > 0:  # Only label bars with non-zero values
            plt.text(i + 1, v, str(v), ha='center', va='bottom', fontsize=8)  # Reduced font size
    
    plt.tight_layout()
    plt.savefig('haplotype_1_frequency_distribution_.png')
    # plt.savefig('haplotype_2_frequency_distribution_.png')
    plt.close()

    print("Total haplotypes:", sum(y))
    print("Distribution:")
    for count, freq in sorted(counter.items()):
        print(f"{count} sample(s): {freq} haplotype(s)")

# Main execution
file_path = '/gs/gsfs0/shared-lab/greally-lab/David/simple_allele-stacker/filtered_consensus_regions/filtered_consensus_H1.bed'
#file_path = './filtered_consensus_H2_cleaned.bed'
sample_counts = parse_bed_file(file_path)
plot_frequency_distribution(sample_counts)