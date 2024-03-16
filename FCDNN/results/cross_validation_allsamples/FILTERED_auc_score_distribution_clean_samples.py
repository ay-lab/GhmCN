# Sconda conda activate abc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/; python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Load the CSV file into a DataFrame
file_path = './FILTERED_cross_validation_all_scores_per_sample_per_chrom.csv'
auc_data = pd.read_csv(file_path)

# We need to aggregate AUC values by sample across all chromosomes.
# Group the data by sample and then list all AUC values
sample_auc_data = auc_data.groupby('sample')['auc'].apply(list).reset_index(name='auc_values')

# Explode the auc_values so each value has its own row (required for Seaborn's boxplot)
sample_auc_data = sample_auc_data.explode('auc_values')

# Convert auc_values back to float (explode makes them objects)
sample_auc_data['auc_values'] = sample_auc_data['auc_values'].astype(float)

# Create the boxplot with Seaborn
plt.figure(figsize=(14, 8))
sns.set(style="whitegrid")
ax = sns.boxplot(x='sample', y='auc_values', data=sample_auc_data, color='lightblue')
sns.stripplot(x='sample', y='auc_values', data=sample_auc_data, color='black', size=4, jitter=True, alpha=0.6)


# Customizing the plot for publication
ax.set_title('AUC Distribution per Sample Across Chromosomes', fontsize=24, fontname='Times New Roman')
ax.set_xlabel('Sample', fontsize=20, fontname='Times New Roman')
ax.set_ylabel('AUC Value', fontsize=20, fontname='Times New Roman')
ax.tick_params(axis='x', labelsize=10, rotation=90)
ax.tick_params(axis='y', labelsize=14)

# Tight layout for better spacing
plt.tight_layout()

# Save the figure
output_sample_boxplot_path = './FILTERED_cross_validation_auc_distribution_per_sample_boxplot'
plt.savefig(output_sample_boxplot_path+'.png', dpi=300, bbox_inches='tight')
plt.savefig(output_sample_boxplot_path+'.pdf', dpi=150, bbox_inches='tight')

# Show the plot
plt.show()
plt.close()












# Add red lines
# Merging the new AUC values into the original AUC distribution data
# This merge will allow us to know the specific AUC value to draw the horizontal line for each sample
# Load the new CSV file containing the specific AUC values per sample
new_auc_file_path = './FILTERED_auc_score_distribution_clean_samples.csv'
new_auc_data = pd.read_csv(new_auc_file_path)
merged_auc_data = auc_data.merge(new_auc_data[['sample', 'auc']], on='sample', suffixes=('', '_specific'))
model_type = new_auc_data['Model'].iloc[0]  # Assuming all rows have the same model type for simplicity



# Now let's create the boxplot with a red horizontal line for each sample
plt.figure(figsize=(16, 10))
sns.set(style="whitegrid")
# Creating the boxplot
sns.stripplot(x='sample', y='auc', data=merged_auc_data, color='black', size=4, jitter=True, alpha=0.6)
ax = sns.boxplot(x='sample', y='auc', data=merged_auc_data, color='lightblue')
# Adding a red horizontal line for each sample's specific AUC value
unique_samples = merged_auc_data['sample'].unique()
# To correctly center the red lines over each boxplot, we'll calculate the positions more accurately
for i, sample in enumerate(unique_samples):
    specific_auc = merged_auc_data.loc[merged_auc_data['sample'] == sample, 'auc_specific'].values[0]
    # The position of the sample's boxplot is determined by its order in the plot
    plt.plot([i - 0.4, i + 0.4], [specific_auc, specific_auc], color='red', linestyle='-', linewidth=2)

# To add a correctly colored red line in the legend, we use a workaround by plotting an invisible line
# and then creating a legend for it
_, = plt.plot([], [], color='red', label=f'{model_type} Model AUC value ')
plt.legend(loc='lower left', fontsize=14, frameon=True, handlelength=2, handletextpad=0.5, edgecolor="black", fancybox=True)

# Customizing the plot
ax.set_title('AUC Distribution per Sample Across Chromosomes with Specific AUC', fontsize=20)
ax.set_xlabel('Sample', fontsize=16)
ax.set_ylabel('AUC Value', fontsize=16)
ax.tick_params(axis='x', labelsize=10, rotation=90)
ax.tick_params(axis='y', labelsize=14)

plt.tight_layout()

# Save the figure
output_sample_boxplot_path = './FILTERED_cross_validation_auc_distribution_per_sample_boxplot_overlayed'
plt.savefig(output_sample_boxplot_path+'.png', dpi=300, bbox_inches='tight')
plt.savefig(output_sample_boxplot_path+'.pdf', dpi=150, bbox_inches='tight')

# Show the plot
plt.close()











# Per chromosome
# Plotting AUC distribution chromosome-wise with a red line representing the mean or median of the new data AUC scores
# First, let's decide whether to use the mean or median for the red line. For this example, we'll use the mean.
import re

mean_auc_score = 0.87
# Re-defining the sort function with the 're' module imported
def sort_chromosomes(chromosomes):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(chromosomes, key=alphanum_key)

# Applying the custom sort function to the chromosome names again
sorted_chromosomes = sort_chromosomes(auc_data['chr'].unique())

# Now, creating the boxplot chromosome-wise
plt.figure(figsize=(16, 10))
sns.set(style="whitegrid")

# Creating the boxplot with chromosomes sorted correctly
sns.stripplot(x='chr', y='auc', data=auc_data, color='black', order=sorted_chromosomes, size=4, jitter=True, alpha=0.6)
ax = sns.boxplot(x='chr', y='auc', data=auc_data, color='lightblue', order=sorted_chromosomes)

# Adding black points for each data point on top of the boxplots

# Adding a red line for the mean AUC score
plt.axhline(mean_auc_score, color='red', linestyle='--', linewidth=1, label=f'{model_type} Model mean AUC: {mean_auc_score:.2f}')

# Adding a legend
plt.legend(loc='upper right', fontsize=14, frameon=True)

# Customizing the plot
ax.set_title('AUC Distribution per Chromosomes', fontsize=20)
ax.set_xlabel('Chromosome', fontsize=16)
ax.set_ylabel('AUC Value', fontsize=16)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)

plt.tight_layout()

# Save the figure with the mean AUC line
output_chrom_wise_with_mean_path = './FILTERED_cross_validation_auc_distribution_per_chrom_boxplot_overlayed'
plt.savefig(output_chrom_wise_with_mean_path+'.png', dpi=300, bbox_inches='tight')
plt.savefig(output_chrom_wise_with_mean_path+'.pdf', dpi=150, bbox_inches='tight')

# Show the plot
plt.show()
plt.close()
