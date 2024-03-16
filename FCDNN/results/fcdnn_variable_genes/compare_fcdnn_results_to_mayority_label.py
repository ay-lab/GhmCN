# Sconda conda activate abc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_variable_genes; python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc
from scipy.stats import shapiro, mannwhitneyu
import os
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
acc_columns = ['gene', 'tpm']
auc_columns = ['gene', 'category']
for i, name in enumerate(samples.tolist()):
    dev_m = pd.read_table(f'{path_in}/{name}_Dev_M.txt', names=colnames)[['gene', 'category', 'tpm']]
    tst_m = pd.read_table(f'{path_in}/{name}_Test_M.txt', names=colnames)[['gene', 'category', 'tpm']]
    trn_m = pd.read_table(f'{path_in}/{name}_Train_M.txt', names=colnames)[['gene', 'category', 'tpm']]
    all = pd.concat([dev_m, tst_m, trn_m])
    all = all[~all['gene'].isin(noisy_genes)]
    all_auc = all[auc_columns].copy()
    all = all[acc_columns].copy()
    all_classes = all if i == 0 else pd.merge(all_classes, all, on=["gene"], how="left")
    all_classes_auc = all_auc if i == 0 else pd.merge(all_classes_auc, all_auc, on=["gene"], how="left")
    print(name)

all_classes.columns = ['gene'] + samples.tolist()
all_classes_auc.columns = ['gene'] + samples.tolist()

threshold = all_classes.sum(axis=1).median()/(all_classes.shape[1]-1)
all_genes = all_classes.mean(axis=1)
labels = (all_genes > threshold).astype(int)

# __IMPORTANT__
# For this section below run first 
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_variable_genes/fcdnn_variable_genes.sh
# __^^^^^^^^^__
variable_genes = pd.read_csv('/tmp/all_any_variable_genes.txt', names = ['gene'])
all_classes.reset_index(inplace=True)
all_classes_auc.reset_index(inplace=True)
# Get the mayority vote labes of most variable genes:
variable_genes_indexes = all_classes.reset_index().merge(variable_genes,how='inner', on=['gene'])['index']
labels = labels[all_classes['index'].isin(variable_genes_indexes)]
all_classes = all_classes[all_classes['index'].isin(variable_genes_indexes)].drop(columns=['index'])
all_classes_auc = all_classes_auc[all_classes_auc['index'].isin(variable_genes_indexes)].drop(columns=['index'])
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
    # Get accuracy in excluded sample
    current_sample_acc = current_sample_acc.query('most_frequent_label == category').shape[0] / current_sample_acc.shape[0]
    accuracies_per_sample.append(current_sample_acc)
    # Get AUC in excluded sample
    excluded_auc = all_classes_auc.loc[:, excluded_idx]
    included_auc = all_classes_auc.loc[:, included_idx]
    prob = ((included_auc == "High").sum(axis=1) / (samples.shape[0] - 1)).tolist()
    obs = (excluded_auc.iloc[:,0] == "High").astype(int)
    fpr, tpr, _ = roc_curve(obs, prob)
    auc_per_sample.append(auc(fpr, tpr))
    print(i)


plt.figure(figsize=(6, 10))
sns.boxplot(y=accuracies_per_sample, color='lightblue')
sns.swarmplot(y=accuracies_per_sample, color='darkblue')
plt.title('[Baseline - Any Variable Genes]\nAccuracy of Label Prediction per Sample - TPM', fontsize=16)
plt.xlabel('Samples', fontsize=14)
plt.ylabel('Accuracy', fontsize=14)
plt.ylim(0.5,0.8)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('./FILTERED_baseline_acurracy_per_sample_balanced_split_any_variable_genes.pdf')
plt.close()


plt.figure(figsize=(6, 10))
sns.boxplot(y=auc_per_sample, color='lightblue')
sns.swarmplot(y=auc_per_sample, color='darkblue')
plt.title('[Baseline - Any Variable Genes]\nAUC Distribution. 1 Dot per Sample')
plt.xlabel('Samples')
plt.ylim(0.7,1.0)
plt.ylabel('AUC')
plt.savefig('FILTERED_baseline_tpm_auc_per_sample_balanced_split_any_variable_genes.pdf')
plt.close()

combined_model = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_variable_genes/predictions_CleanSamples_tested_on_AnyVariableGenes.csv')
combined_model.Predicted = combined_model.Predicted.astype(int)
samples=49
genes_per_sample =  combined_model.shape[0]//samples

combined_acc = []
for i in range(samples):
    chunk = combined_model[i*genes_per_sample:(i+1)*genes_per_sample].copy()
    combined_acc.append(chunk.query('Observed  == Predicted').shape[0] / chunk.shape[0])

merged_data = pd.DataFrame({'Model':combined_acc, 'Baseline':accuracies_per_sample })
n_model = len(merged_data['Model'])
n_baseline = len(merged_data['Baseline'])

max_values = {0:merged_data['Model'].max(),1:merged_data['Baseline'].max()}

# Shapiro-Wilk test, which is what `shapiro` function performs, calculates a test statistic
# that represents likelyhood of a given dataset being drawn from a Gaussian distribution.
stat1, p1 = shapiro(merged_data['Model'])
stat2, p2 = shapiro(merged_data['Baseline'])
# Decide which test to use based on the p-value
if p1 < 0.05 or p2 < 0.05:
    print("At least one of the groups is not normally distributed. Use Mann-Whitney U test.\n")
    print(f'P values:\n\tModel {round(p1,4)}\n\tBaseline {round(p2,4)}')
else:
    print("Both groups are normally distributed. You may consider using t-test.")

# Assuming `group1` and `group2` are the data for the two groups
stat, p_value = mannwhitneyu(merged_data['Model'], merged_data['Baseline'], alternative='two-sided')
print(f'Statistic: {stat}, P-value: {p_value}')




# Merged boxplot:
font_name = 'Nimbus Roman'
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=merged_data, width=0.5, palette="Blues", fliersize=5)
sns.swarmplot(data=merged_data, size=4, edgecolor="gray", linewidth=0.5)
ax.set_title('Baseline and Combined Model Accuracy of Label Prediction\nIn Genes With Any Expression Change', fontsize=24, fontname=font_name)
ax.set_xlabel('', fontsize=24, fontname=font_name)
ax.set_ylabel('Accuracy', fontsize=24, fontname=font_name)
plt.text(0, float(max_values[0]), f'N={n_model}', ha='center', va='bottom',fontname=font_name, fontsize=14)
plt.text(1, float(max_values[1]), f'N={n_baseline}', ha='center', va='bottom',fontname=font_name, fontsize=14)
plt.xticks(fontsize=22, fontname=font_name)
plt.yticks(fontsize=22, fontname=font_name)
plt.yticks(np.arange(0.65, 0.85+0.05, 0.05))
plt.savefig('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_variable_genes/FILTERED_Baseline_and_Combined_Model_Accuracy_of_Label_Prediction.pdf')
plt.close()