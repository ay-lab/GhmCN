# Sconda conda activate abc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/genes_coef_var_sq; python
# python expression_threshold_per_celltype.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


# Function to TPM normalize the gene counts
def tpm_normalize(df):
    # Assuming column 6 has gene length in base pairs (bp)
    length_kb = df.iloc[:, 5] / 1000  # Convert gene length from bp to kilobases (kb)
    # Assuming column 7 has raw read counts
    reads = df.iloc[:, 6]
    # Calculate Reads Per Kilobase (RPK)
    rpk = reads / length_kb
    # Calculate the scaling factor by summing all RPK values and dividing by 1,000,000
    scaling_factor = rpk.sum() / 1_000_000
    # Divide RPK by the scaling factor to get TPM
    tpm = rpk / scaling_factor
    return tpm

base_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/gene_expression_class_concordance'
out_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/gene_expression_distribution'
samples = pd.read_csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/filtered_files_surnames.csv", header = None, usecols = [7,10])
samples.columns = ["sample_name","rna_file"]

all_normalized_expression = []
for file in samples.rna_file:
    feat_counts = pd.read_table(file, comment="#")
    tpm_values = tpm_normalize(feat_counts)
    all_normalized_expression.append(tpm_values)

all_normalized_expression = pd.concat(all_normalized_expression, axis=1, ignore_index=True)
all_normalized_expression.index = feat_counts['Geneid'].values

mean_tpm = all_normalized_expression.mean(axis=1)
std_tpm = all_normalized_expression.std(axis=1)
cv2 = (std_tpm / mean_tpm) ** 2

# Plot the log10(CV2) in the Y axis versus log10(mean TPM expression) in X
plt.figure(figsize=(10, 6))
sns.scatterplot(x=np.log10(mean_tpm), y=np.log10(cv2), alpha=0.5)
# Add a loess regression lines fitting the dispersion to the mean expression level
sns.regplot(x=np.log10(mean_tpm), y=np.log10(cv2), scatter=False, lowess=True, color='red')

plt.xlabel('mean expression [log10(TPM)]')
plt.ylabel('variability [log10(CV²)]')
# plt.title('log10(CV²) vs log10(Mean TPM) per Gene')
plt.title('Variability across gene expression TPM level')
plt.savefig('FILTERED_gene_cv2_with_loess.png')
plt.savefig('FILTERED_gene_cv2_with_loess.pdf')

