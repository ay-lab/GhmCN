# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/genes_per_category; python
# Check total genes per category for the three classes:
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
C=2
split = 'Balanced' # 'Practical'
data_path = os.path.join('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/10.AutoFormattedDatasets/02.Datasets',f'Classes{C}_'+split)
colnames = ['gene','catergory','digit_category','chrom','index','tpm']
samples = pd.read_csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/filtered_files_surnames.csv", header = None)[7]
samples.columns = ["sample_name","rna_file"]

n_off = []
n_off25p = []
n_off33p = []
n_on_25p = []
n_on_33p = []
n_low_25p = []
n_high_25p = []
n_low_50p = []
n_high_50p = []
for sample_name in samples.tolist():
    dev_m = pd.read_table(data_path+f'/{sample_name}_Dev_M.txt',names = colnames)
    tst_m = pd.read_table(data_path+f'/{sample_name}_Test_M.txt',names = colnames)
    trn_m = pd.read_table(data_path+f'/{sample_name}_Train_M.txt',names = colnames)
    all = pd.concat([dev_m, tst_m, trn_m])
    expressed = all.query('tpm > 0')
    threshold_25p = np.percentile(expressed.tpm,25)
    threshold_33p = np.percentile(expressed.tpm,33)
    threshold_50p = np.percentile(expressed.tpm,50)
    # off
    n_off.append(all.query('tpm == 0').shape[0])
    n_off25p.append(all.query('tpm == 0').shape[0] + sum(expressed.tpm <= threshold_25p))
    n_off33p.append(all.query('tpm == 0').shape[0] + sum(expressed.tpm <= threshold_33p))
    # on
    n_on_25p.append(sum(expressed.tpm > threshold_25p))
    n_on_33p.append(sum(expressed.tpm > threshold_33p))
    # low
    n_low_25p.append(sum(expressed.tpm <= threshold_25p))
    n_low_50p.append(sum(expressed.tpm <= threshold_50p))
    # high
    n_high_25p.append(sum(expressed.tpm > threshold_25p))
    n_high_50p.append(sum(expressed.tpm > threshold_50p))

df2_25p = pd.DataFrame({ 'off' : n_off25p, 'on' :  n_on_25p})
df2_33p = pd.DataFrame({ 'off' : n_off33p, 'on' :  n_on_33p})

df_25p = pd.DataFrame({ 'off' : n_off, 'low': n_low_25p, 'high' :  n_high_25p})
df_50p = pd.DataFrame({ 'off' : n_off, 'low': n_low_50p, 'high' :  n_high_50p})

plt.boxplot([df_25p['off'], df_25p['low'], df_25p['high']])
plt.xticks([1, 2, 3], ['off', 'low', 'high'])
plt.xlabel('Category')
plt.ylabel('# genes')
plt.title('Total Genes per Category @ 25% Threshold')
plt.savefig('./FILTERED_genes_per_category_25p_ternary.png')
plt.close()

plt.boxplot([df_50p['off'], df_50p['low'], df_50p['high']])
plt.xticks([1, 2, 3], ['off', 'low', 'high'])
plt.xlabel('Category')
plt.ylabel('# genes')
plt.title('Total Genes per Category @ 50% Threshold')
plt.savefig('./FILTERED_genes_per_category_50p_ternary.png')
plt.close()

plt.boxplot([df2_25p['off'], df2_25p['on']])
plt.xticks([1, 2], ['off', 'on'])
plt.xlabel('Category')
plt.ylabel('# genes')
plt.title('Total Genes per Category @ 25% Threshold - Binary')
plt.savefig('./FILTERED_genes_per_category_25p_binary.png')
plt.close()

plt.boxplot([df2_33p['off'], df2_33p['on']])
plt.xticks([1, 2], ['off', 'on'])
plt.xlabel('Category')
plt.ylabel('# genes')
plt.title('Total Genes per Category @ 33% Threshold - Binary')
plt.savefig('./FILTERED_genes_per_category_33p_binary.png')
plt.close()