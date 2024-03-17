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
    dev_m = pd.read_table(f'{path_in}/{name}_Dev_M.txt', names=colnames)[['gene', 'tpm']]
    tst_m = pd.read_table(f'{path_in}/{name}_Test_M.txt', names=colnames)[['gene', 'tpm']]
    trn_m = pd.read_table(f'{path_in}/{name}_Train_M.txt', names=colnames)[['gene', 'tpm']]
    all = pd.concat([dev_m, tst_m, trn_m])
    all = all[~all['gene'].isin(noisy_genes)]
    all_classes = all if i == 0 else pd.merge(all_classes, all, on=["gene"], how="left")
    print(name)

all_classes.columns = ['gene'] + samples.tolist()
threshold = all_classes.sum(axis=1).median()/(all_classes.shape[1]-1)
all_genes = all_classes.mean(axis=1) # this is a neat trick because not only calculates the mean, it repeats its value across all rows!
labels = (all_genes > threshold).astype(int) # That's why I can use this operation here, each gene (column-wise) has the mean, and is 

all_classes_bkup = all_classes.copy()
all_classes = all_classes_bkup.copy()

accuracies_per_sample = []
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
    sample_threshold = excluded.median(axis = 0)[0]
    excluded = (excluded > sample_threshold).astype(int)
    
    
    # Check for excluded sample names:
    # [elem for elem, include in zip(samples.tolist(), excluded_idx) if include]
    
    # from excluded data type keep only consistent genes:
    # excluded = excluded.loc[excluded.drop('gene', axis=1).nunique(axis=1) == 1].iloc[:,[0,1]]
    # excluded.columns = ['gene', "category"]
    
    # also keep only consistent from all the other datasts:
    # included = included[included.gene.isin(excluded.gene.tolist())]
    
    # Get most frequent label:
    # most_frequent_label = included.mode(axis=1)[0]
    current_sample_acc = pd.DataFrame({
        'most_frequent_label': labels,
        'category': excluded.iloc[:,0]
    })
    fpr, tpr, _ = roc_curve(excluded.iloc[:,0], labels)
    auc_per_sample.append(auc(fpr, tpr))
    # Get accuracy in excluded sample
    current_sample_acc = current_sample_acc.query('most_frequent_label == category').shape[0] / current_sample_acc.shape[0]
    accuracies_per_sample.append(current_sample_acc)
    print(i)


plt.figure(figsize=(6, 10))
sns.boxplot(y=accuracies_per_sample, color='lightblue')
sns.swarmplot(y=accuracies_per_sample, color='darkblue')
plt.title('[Baseline] Accuracy of Label Prediction per Sample')
plt.xlabel('Samples')
plt.ylabel('Accuracy')
plt.ylim(0.7,1.0)
plt.savefig('FILTERED_baseline_tpm_acurracy_per_sample_balanced_split.pdf')
plt.close()

# new_auc = pd.read_csv('clean_data_auc.txt', header = None)[0]
# plt.figure(figsize=(6, 10))
# sns.boxplot(y=new_auc, color='lightblue')
# sns.swarmplot(y=new_auc, color='darkblue')
# plt.title('AUC Score of Label Prediction per Sample')
# plt.xlabel('Samples')
# plt.ylabel('AUC')
# plt.ylim(0.7,1.0)
# plt.savefig('FILTERED_auc_per_sample_balanced_split.pdf')
# plt.close()

# new_accuracy = pd.read_csv('clean_data_accuracy.txt', header = None)[0]
# plt.figure(figsize=(6, 10))
# sns.boxplot(y=new_accuracy, color='lightblue')
# sns.swarmplot(y=new_accuracy, color='darkblue')
# plt.title('Accuracy of Label Prediction per Sample')
# plt.xlabel('Samples')
# plt.ylabel('AUC')
# plt.ylim(0.7,1.0)
# plt.savefig('FILTERED_accuracy_per_sample_balanced_split.pdf')
# plt.close()


# For this section below run first 
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_most_variable_genes/fcdnn_most_variable_genes.sh

variable_genes = pd.read_csv('/tmp/all.txt', names = ['gene'])
all_classes.reset_index(inplace=True)
all_classes.shape
# Get the mayority vote labes of most variable genes:
labels = labels[all_classes['index'].isin(all_classes.reset_index().merge(variable_genes,how='inner', on=['gene'])['index'])]
all_classes = all_classes[all_classes['index'].isin(all_classes.reset_index().merge(variable_genes,how='inner', on=['gene'])['index'])].drop(columns=['index'])
all_classes.shape
# all_classes = all_classes_bkup.copy()
# all_classes_bkup = all_classes.copy()

accuracies_per_sample = []
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
    sample_threshold = excluded.median(axis = 0)[0]
    excluded = (excluded > sample_threshold).astype(int)
    
    # Get most frequent label:
    # most_frequent_label = included.mode(axis=1)[0]
    current_sample_acc = pd.DataFrame({
        'most_frequent_label': labels,
        'category': excluded.iloc[:,0]
    })
    fpr, tpr, _ = roc_curve(excluded.iloc[:,0], labels)
    auc_per_sample.append(auc(fpr, tpr))
    # Get accuracy in excluded sample
    current_sample_acc = current_sample_acc.query('most_frequent_label == category').shape[0] / current_sample_acc.shape[0]
    accuracies_per_sample.append(current_sample_acc)
    print(i)


plt.figure(figsize=(6, 10))
sns.boxplot(y=accuracies_per_sample, color='lightblue')
sns.swarmplot(y=accuracies_per_sample, color='darkblue')
plt.title('[Baseline] Accuracy of Label Prediction per Sample\nMost Variable Genes (2446)')
plt.xlabel('Samples')
plt.ylabel('Accuracy')
plt.ylim(0.4,1.0)
plt.savefig('FILTERED_baseline_tmp_acurracy_per_sample_balanced_split_most_variable_genes_2446.pdf')
plt.close()

# new_auc = pd.read_csv('clean_data_auc.txt', header = None)[0]
plt.figure(figsize=(6, 10))
sns.boxplot(y=auc_per_sample, color='lightblue')
sns.swarmplot(y=auc_per_sample, color='darkblue')
plt.title('[Baseline] AUC of Label Prediction per Sample\nMost Variable Genes (2446)')
plt.xlabel('Samples')
plt.ylabel('AUC')
plt.ylim(0.0,1.0)
plt.savefig('FILTERED_baseline_tmp_auc_per_sample_balanced_split_most_variable_genes_2446.pdf')
plt.close()

