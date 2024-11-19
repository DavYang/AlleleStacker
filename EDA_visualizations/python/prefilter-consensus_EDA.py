import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import ScalarFormatter, FuncFormatter

# Read the TSV file
df = pd.read_csv('consensus_regions_region_statistics.tsv', sep='\t')

# Set up the plot style
plt.style.use('ggplot')
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 25))

# Function to format y-axis labels in scientific notation
def scientific(x, pos):
    return f'${x:.0e}$'

# 1. Bar plot of H1 and H2 counts without error bars
sns.barplot(x='Sample', y='value', hue='variable', data=pd.melt(df, id_vars=['Sample'], value_vars=['H1 Count', 'H2 Count']), 
            palette={'H1 Count': 'b', 'H2 Count': 'r'}, alpha=0.5, ax=ax1)

ax1.set_title('H1 and H2 Counts by Sample')
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
ax1.set_xlabel('Sample')
ax1.set_ylabel('Count')
ax1.yaxis.set_major_formatter(FuncFormatter(scientific))
ax1.legend(title='Haplotype')

# 2. Scatter plot of H1 Mean Size vs H2 Mean Size
sns.scatterplot(x='H1 Mean Size', y='H2 Mean Size', data=df, ax=ax2)
ax2.set_title('H1 Mean Size vs H2 Mean Size')
ax2.set_xlabel('H1 Mean Size (base pairs)')
ax2.set_ylabel('H2 Mean Size (base pairs)')
ax2.xaxis.set_major_formatter(FuncFormatter(scientific))
ax2.yaxis.set_major_formatter(FuncFormatter(scientific))
for i, sample in enumerate(df['Sample']):
    ax2.annotate(sample, (df['H1 Mean Size'][i], df['H2 Mean Size'][i]))

# 3. Box plot of H1 and H2 Size Variances
df_variance = df.melt(id_vars=['Sample'], value_vars=['H1 Size Variance', 'H2 Size Variance'], var_name='Haplotype', value_name='Size Variance')
sns.boxplot(x='Haplotype', y='Size Variance', data=df_variance, ax=ax3)
ax3.set_title('Distribution of H1 and H2 Size Variances')
ax3.set_xlabel('Haplotype')
ax3.set_ylabel('Size Variance (base pairs²)')
ax3.yaxis.set_major_formatter(FuncFormatter(scientific))

# 4. New plot: Bar plot of Mean Region Size with Variance as error bars
df_mean_var = df.melt(id_vars=['Sample'], 
                      value_vars=['H1 Mean Size', 'H2 Mean Size', 'H1 Size Variance', 'H2 Size Variance'],
                      var_name='Metric', value_name='Value')
df_mean_var['Haplotype'] = df_mean_var['Metric'].str[:2]
df_mean_var['Type'] = df_mean_var['Metric'].str[3:]

df_mean = df_mean_var[df_mean_var['Type'] == 'Mean Size'].copy()
df_var = df_mean_var[df_mean_var['Type'] == 'Size Variance'].copy()

df_mean['Variance'] = df_var['Value'].values

sns.barplot(x='Sample', y='Value', hue='Haplotype', data=df_mean, 
            palette={'H1': 'b', 'H2': 'r'}, alpha=0.5, ax=ax4)

# Add error bars
for i, sample in enumerate(df['Sample']):
    for j, haplotype in enumerate(['H1', 'H2']):
        data = df_mean[(df_mean['Sample'] == sample) & (df_mean['Haplotype'] == haplotype)]
        ax4.errorbar(x=i + (j-0.5)*0.4, y=data['Value'], yerr=np.sqrt(data['Variance']),
                     fmt='none', c='k', capsize=3, alpha=0.7)

ax4.set_title('Mean Region Size by Sample and Haplotype (with Size Variance)')
ax4.set_xticklabels(ax4.get_xticklabels(), rotation=45, ha='right')
ax4.set_xlabel('Sample')
ax4.set_ylabel('Mean Region Size (base pairs)')
ax4.yaxis.set_major_formatter(FuncFormatter(scientific))
ax4.legend(title='Haplotype')

plt.tight_layout()
plt.savefig('consensus_regions_visualization.png', dpi=300, bbox_inches='tight')
plt.show()