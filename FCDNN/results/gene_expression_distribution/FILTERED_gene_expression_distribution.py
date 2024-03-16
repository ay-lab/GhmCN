# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/gene_expression_distribution; python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Calculate the mean and standard deviation of the TPM values
def plot_tmp_dist(tpm_values:list, name:str, split_thresholds:pd.DataFrame, is_biol:bool = False, is_log:bool = True):
    scale = 'log10_tpm' if is_log else 'raw_tpm'
    split = 'Meaningful' if is_biol else 'Balanced'
    mean_tpm = np.mean(tpm_values, axis=0)
    std_tpm = np.std(tpm_values, axis=0)
    # median thresholds:
    split_idx = 0 if is_biol else 1
    median_threshold = split_thresholds.iloc[split_idx,1]
    bottom_q = split_thresholds.iloc[split_idx,2] # 2.5%
    top_q = split_thresholds.iloc[split_idx,5] # 97.5%
    if is_log:
        median_threshold, bottom_q, top_q = np.log10(median_threshold), np.log10(bottom_q), np.log10(top_q)
    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.plot(mean_tpm, label='Mean Expression')
    plt.fill_between(range(len(mean_tpm)), mean_tpm - std_tpm, mean_tpm + std_tpm, alpha=0.2, label='Standard Deviation')
    plt.xlabel('Expression Rank')
    plt.ylabel('Expression Log10(TPM)')
    plt.title(f'Sorted Mean Expression Distribution with Standard Deviation - {name}')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=1)
    plt.axhline(y=median_threshold, color='orange', linestyle='--', linewidth=1, label = f'Median Threshold @ {split} Data Split')
    plt.fill_between(range(len(mean_tpm)), bottom_q, top_q, alpha=0.2, color='orange', label='+/-95% Thresholds')
    plt.legend()
    plt.savefig(f'./FILTERED_gene_expression_dist_{name}_{scale}_{split}.pdf')


def plot_mean_var(tpm_values:list, name:str, is_log:bool = True):
    scale = 'log10_tpm' if is_log else 'raw_tpm'
    scale_label = ' (Log10 TPM)' if is_log else ' (Raw TPM)'
    mean_tpm = np.mean(tpm_values, axis=0)
    variance_log10_tpm = np.var(tpm_values, axis=0)
    # Sort the data by mean log10(TPM) values
    # sorted_data = sorted(zip(mean_tpm, variance_log10_tpm), key=lambda x: x[0])
    # sorted_mean_log10_tpm, sorted_variance_log10_tpm = zip(*sorted_data)
    # Create a line plot for the mean-variance relationship
    plt.figure(figsize=(10, 6))
    plt.plot(mean_tpm, variance_log10_tpm, label='Variance')
    plt.xlabel('Expression Rank Mean Expression' + scale_label)
    plt.ylabel('Variance of Log10 Expression'+scale_label)
    plt.title(f'Sorted Mean versus Variance - {name}')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=1)
    plt.legend()
    plt.savefig(f'./FILTERED_gene_expression_mean_var_{name}_{scale}.pdf')


# Load the median threshold per split
out_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/gene_expression_distribution'
out_file_overall=os.path.join(out_path,'FILTERED_median_threshold_per_split.csv')
split_thresholds = pd.read_csv(out_file_overall)


C=2
split = 'Balanced' # 'Practical'
fcdnn='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN'
data_path = os.path.join('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets',f'Classes{C}_'+split)
colnames = ['gene','catergory','digit_category','chrom','index','tpm']
all_data = pd.read_csv(f'{fcdnn}/filtered_files_surnames.csv', header = None)
samples = pd.read_csv(f'{fcdnn}/filtered_files_surnames.csv', header = None)[7]
sample_groups = pd.read_csv(f'{fcdnn}/data/sample_groups.csv', header = None, names = ['sample','group'])

# Create an empty list to store the 'tpm' values from each file
tpm_values = []
for name in samples.tolist():
    dev_m = pd.read_table(data_path+f'/{name}_Dev_M.txt',names = colnames)
    tst_m = pd.read_table(data_path+f'/{name}_Test_M.txt',names = colnames)
    trn_m = pd.read_table(data_path+f'/{name}_Train_M.txt',names = colnames)
    # Extract the 'tpm' column, sort values and append it to tpm_values
    tpm_values.append(np.log10(pd.concat([dev_m, tst_m, trn_m])['tpm'].sort_values() + 1e-10))


# plot_mean_var(tpm_values = tpm_values, name = 'All_Samples', is_log = True)
plot_tmp_dist(tpm_values = tpm_values, name = 'Clean_Samples', split_thresholds=split_thresholds, is_biol = False, is_log = True)
plot_tmp_dist(tpm_values = tpm_values, name = 'Clean_Samples', split_thresholds=split_thresholds, is_biol = True, is_log = True)
plot_mean_var(tpm_values = tpm_values, name = 'Clean_Samples', is_log = True)
plot_mean_var(tpm_values = tpm_values, name = 'Clean_Samples', is_log = False)
