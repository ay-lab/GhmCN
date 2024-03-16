# Take 1 cell-type, say inkt
# Take its metadata for everything: training, testing, dev
# Make a pandas df with gene and category

# Then for each of the other samples, join gene and category.
# merge by gene

# Sconda conda activate TensorFlow_CPU_01
# cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/baseline_compare_mayority_of_labels; python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.metrics import roc_curve, auc

C=2
base_path = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data'
path_in = os.path.join(base_path,f'C{C}B')
colnames = ['gene','category','digit_category','chrom','index','tpm']
samples = pd.read_csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/filtered_files_surnames.csv", header = None)[7]
noisy_genes = [
    "Mir5098",
    "Eno1b",
    "Mir684-1",
    "Gm5643",
    "Gm5801",
    "Bc1",
    "Gm5512",
    "Btg3",
    "Mir3473a",
]

# Improved version of the provided code

# Load all data
for i, name in enumerate(samples.tolist()):
    dev_m = pd.read_table(f'{path_in}/{name}_Dev_M.txt', names=colnames)[['gene', 'category']]
    tst_m = pd.read_table(f'{path_in}/{name}_Test_M.txt', names=colnames)[['gene', 'category']]
    trn_m = pd.read_table(f'{path_in}/{name}_Train_M.txt', names=colnames)[['gene', 'category']]
    all = pd.concat([dev_m, tst_m, trn_m])
    all = all[~all['gene'].isin(noisy_genes)]
    all_classes = all if i == 0 else pd.merge(all_classes, all, on=["gene"], how="left")
    print(name)

all_classes.columns = ['gene'] + samples.tolist()

auc_per_sample = []
excl_idx = [False] + [False] * (len(samples))
incl_idx = [False] + [True] * (len(samples))
for i, sample_out in enumerate(samples.tolist()):
    excluded_idx = [False] + [False] * (len(samples))
    included_idx = [False] + [True] * (len(samples))
    excluded_idx[i+1] = True
    included_idx[i+1] = False
    excluded = all_classes.loc[:, excluded_idx]
    included = all_classes.loc[:, included_idx]
    
    # Check for excluded sample names:
    # [elem for elem, include in zip(samples.tolist(), excluded_idx) if include]
    
    # from excluded data type keep only consistent genes:
    # excluded = excluded.loc[excluded.drop('gene', axis=1).nunique(axis=1) == 1].iloc[:,[0,1]]
    # excluded.columns = ['gene', "category"]
    
    # also keep only consistent from all the other datasts:
    # included = included[included.gene.isin(excluded.gene.tolist())]
    prob = ((included == "High").sum(axis=1) / (samples.shape[0] - 1)).tolist()
    obs = (excluded.iloc[:,0] == "High").astype(int)
    fpr, tpr, _ = roc_curve(obs, prob)
    auc_per_sample.append(auc(fpr, tpr))
    print(i)


plt.figure(figsize=(6, 10))
sns.boxplot(y=auc_per_sample, color='lightblue')
sns.swarmplot(y=auc_per_sample, color='darkblue')
plt.title('[Baseline] AUC Score of Label Prediction per Sample - TPM')
plt.xlabel('Samples')
plt.ylim(0.7,1.0)
plt.ylabel('AUC')
plt.savefig('FILTERED_baseline_tpm_auc_per_sample_balanced_split.pdf')
plt.close()

# Take much smaller datasets by reducing the amount of samples per laboratory/paper
# If multiple celltype, take 1 replicate.
# pd.DataFrame(auc_per_sample).to_csv("baseline_auc_per_samples.txt", sep="\t", index=False,header=False)

auc_per_sample