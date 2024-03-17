# Sconda conda activate abc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/gene_expression_class_concordance
# python rna_replicates_concordance.py
import pandas as pd
import numpy as np
from sklearn.metrics import cohen_kappa_score
import os
import seaborn as sns
import matplotlib.pyplot as plt


# Function to TPM normalize the gene counts
def tpm_normalize(df):
    length = df.iloc[:, 5]  # Assuming column 6 has gene length
    reads = df.iloc[:, 6]   # Assuming column 7 has raw reads
    tpm = (reads / length) * 1e6 / (reads.sum() / 1e6)
    return tpm

# Function to classify genes as lower/higher based on the median threshold
def classify_genes(tpm_values, biologic_split=True):
    non_zero_genes = tpm_values[tpm_values > 0]
    class_threshold = np.median(non_zero_genes) \
        if biologic_split \
        else np.median(tpm_values)
    return np.where(tpm_values <= class_threshold, 'lower', 'higher')

# Function to calculate concordance
def calculate_concordance(labels1, labels2):
    return cohen_kappa_score(labels1, labels2)

# Read the CSV file
base_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/gene_expression_class_concordance'
csv_file = os.path.join(base_path,'rna-reps.txt')
csv_data = pd.read_csv(csv_file)

for biologic_split in [True, False]:
    # Output file
    suffix='biologic_split' if biologic_split else 'balanced_split'
    concordance_file=os.path.join(base_path,f'concordance_results_{suffix}.csv')
    # Iterate through each cell type
    for cell_type in csv_data['cell_type'].unique():
        cell_type_data = csv_data[csv_data['cell_type'] == cell_type]
        # Initialize lists to store labels and concordance values
        all_labels = []
        all_concordance = []
        # Iterate through each replicate
        for _, row in cell_type_data.iterrows():
            # Load and process tab-separated file
            tab_sep_data = pd.read_csv(row['rna_file'], sep='\t', comment='#')
            tpm_values = tpm_normalize(tab_sep_data)
            labels = classify_genes(tpm_values, biologic_split)
            all_labels.append(labels)
        # Calculate concordance across replicates
        for i in range(len(all_labels) - 1):
            for j in range(i + 1, len(all_labels)):
                concordance = calculate_concordance(all_labels[i], all_labels[j])
                all_concordance.append(concordance)
        # Print concordance per cell type
        avg_concordance = np.mean(all_concordance)
        print(f"{cell_type}: {avg_concordance}")
        result_df = pd.DataFrame({'Cell_Type': [cell_type],
                                'Total_Replicates': [len(all_labels)],
                                'Average_Concordance': [np.mean(all_concordance)]})
        result_df.to_csv(concordance_file, mode='a', header=(not os.path.exists(concordance_file)), index=False)
    # Load the concordance results file
    concordance_df = pd.read_csv(concordance_file)
    # Plot and save Boxplot
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=concordance_df, x='Average_Concordance')
    plt.title(f'Boxplot of Average Concordances - {suffix}')
    plt.xlabel('Average Concordance')
    plt.xlim(0.8,1)
    plt.savefig(os.path.join(base_path,f'boxplot_concordances_{suffix}.png'))
    plt.close()
    # Plot and save Violin Plot
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=concordance_df, x='Average_Concordance')
    plt.title(f'Violin Plot of Average Concordances - {suffix}')
    plt.xlabel('Average Concordance')
    plt.xlim(0.75,1.05)
    plt.savefig(os.path.join(base_path,f'violinplot_concordances_{suffix}.png'))
    plt.close()
    # Plot and save Histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(data=concordance_df, x='Average_Concordance', bins=20, kde=True)
    plt.title(f'Histogram of Average Concordances - {suffix}')
    plt.xlabel('Average Concordance')
    plt.xlim(0.8,1)
    plt.savefig(os.path.join(base_path,f'histogram_concordances_{suffix}.png'))
    plt.close()

