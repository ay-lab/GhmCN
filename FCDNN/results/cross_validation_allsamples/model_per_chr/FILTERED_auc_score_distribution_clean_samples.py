# Sconda conda activate abc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/model_per_chr/; python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import json

# Load the CSV file into a DataFrame
out_path = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/model_per_chr'
file_path_chr = os.path.join(out_path,'FILTERED_cross_validation_all_scores_per_chr.tsv')
file_path_sample = os.path.join(out_path,'FILTERED_cross_validation_all_scores_per_sample_per_chrom.tsv')
file_path_genes = os.path.join(out_path,'genes_per_chrom.csv')
# generate a dict kets of the test:dev chromosomes:
dict_file = os.path.join(out_path,'chr_excluded_for_test_and_dev.json')
# Load dict from file
with open(dict_file, 'r') as f:
    dev_chr = json.load(f)

# Handy function to sort data frame based on chromosome:
def sort_df_by_chr(df, dev_chr=dev_chr):
    df['chr_dev'] = df['chr'].map(dev_chr)
    df.dropna(inplace=True)
    df['chr_merged'] = df['chr'].astype(str) + " : " + df['chr_dev'].astype(str)
    df['chr_sort'] = df['chr'].str.extract('(\d+)').astype(int)
    df = df.sort_values('chr_sort').drop(columns=['chr_sort', 'chr_dev', 'chr'])
    return df.copy()


# auc_data_per_chrom = pd.read_csv(file_path)
auc_data = pd.read_csv(file_path_sample)
genes_per_chr = pd.read_csv(file_path_genes)
genes_per_chr = sort_df_by_chr(genes_per_chr)
genes_dict = genes_per_chr.set_index('chr_merged')['gene'].to_dict()

# We need to aggregate AUC values by sample across all chromosomes.
# Group the data by sample and then list all AUC values
sample_auc_data = auc_data.groupby('chr')['auroc'].apply(list).reset_index(name='auc_values')

# Explode the auc_values so each value has its own row (required for Seaborn's boxplot)
sample_auc_data = sample_auc_data.explode('auc_values')
# Sort chromosmes numerically
sample_auc_data = sort_df_by_chr(sample_auc_data)
# Convert auc_values back to float (explode makes them objects)
sample_auc_data['auc_values'] = sample_auc_data['auc_values'].astype(float)

mean_auc_score = 0.869

max_values_by_chr = pd.DataFrame(sample_auc_data.groupby('chr_merged')['auc_values'].max())['auc_values'].to_dict()

# Create the boxplot with Seaborn
plt.figure(figsize=(14, 10))
sns.set_theme(style="whitegrid")
ax = sns.boxplot(x='chr_merged', y='auc_values', data=sample_auc_data, color='lightblue')
sns.stripplot(x='chr_merged', y='auc_values', data=sample_auc_data, color='black', size=4, jitter=True, alpha=0.6)
# Annotating number of genes per chromosome
x_ticks = ax.get_xticks()
chr_labels = [tick.get_text() for tick in ax.get_xticklabels()]
plt.axhline(mean_auc_score, color='red', linestyle='--', linewidth=3, label=f'Combined Model mean AUC: {mean_auc_score:.2f}')
plt.legend(loc='lower left', fontsize=16, frameon=True, handlelength=2, handletextpad=0.5, edgecolor="black", fancybox=True)

# Customizing the plot for publication
font_name = 'Nimbus Roman'
ax.set_title('AUC Distribution per Chromosome Model', fontsize=30, fontname=font_name)
ax.set_xlabel('Chromosomes Excluded From Training\n(Used For Test : Validation)', fontsize=25, fontname=font_name)
ax.set_ylabel('AUC Value', fontsize=25, fontname=font_name)
ax.tick_params(axis='x', labelsize=20, rotation=90)
ax.tick_params(axis='y', labelsize=18)
plt.tight_layout()

# Save the figure
output_sample_boxplot_path = './FILTERED_cross_validation_auc_distribution_per_chr'
plt.savefig(output_sample_boxplot_path+'.pdf', dpi=150, bbox_inches='tight')

# Add genes per chromosome
# Create the boxplot with Seaborn
x_ticks = ax.get_xticks()
chr_labels = [tick.get_text() for tick in ax.get_xticklabels()]
# Iterate through the list of chromosomes to place annotations
for i, chr_label in enumerate(chr_labels):
    gene_count = genes_dict.get(chr_label, 'N/A')  # Fetch gene count for the chromosome
    y_pos = max_values_by_chr.get(chr_label, 1.0)  # Get the upper y-axis limit for the plot
    # Place text annotation slightly above the highest point of the plot
    ax.text(x_ticks[i], y_pos+0.01, f'{gene_count}', ha='center', va='bottom', fontsize=11, color='black',rotation=35, fontname=font_name)

ax.set_ylim(0.63,1.0)
plt.tight_layout()
# Save the figure
output_sample_boxplot_path = './FILTERED_cross_validation_auc_distribution_per_chr_with_n_genes'
plt.savefig(output_sample_boxplot_path+'.pdf', dpi=150, bbox_inches='tight')
plt.close()

