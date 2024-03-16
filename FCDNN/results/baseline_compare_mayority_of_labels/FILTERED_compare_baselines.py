# Take 1 cell-type, say inkt
# Take its metadata for everything: training, testing, dev
# Make a pandas df with gene and category

# Then for each of the other samples, join gene and category.
# merge by gene

# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/baseline_compare_mayority_of_labels; python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import ttest_rel

C=2
base_path = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data'
path_in = os.path.join(base_path,f'C{C}B')
path_in_unbalanced = os.path.join(base_path,f'C{C}P')
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
in_file_class_bal = './all_classes_49_samples_balanced.csv'
in_file_label_bal = './all_classes_49_samples_balanced_labels.csv'

in_file_class_unb = './all_classes_49_samples_unbalanced.csv'
in_file_label_unb = './all_classes_49_samples_unbalanced_labels.csv'

# Load BALANCED data
if not os.path.exists(in_file_class_bal) or not os.path.exists(in_file_label_bal):
    for i, name in enumerate(samples.tolist()):
        dev_m = pd.read_table(f'{path_in}/{name}_Dev_M.txt', names=colnames)[['gene', 'category']]
        tst_m = pd.read_table(f'{path_in}/{name}_Test_M.txt', names=colnames)[['gene', 'category']]
        trn_m = pd.read_table(f'{path_in}/{name}_Train_M.txt', names=colnames)[['gene', 'category']]
        all = pd.concat([dev_m, tst_m, trn_m])
        all = all[~all['gene'].isin(noisy_genes)]
        all_classes = all if i == 0 else pd.merge(all_classes, all, on=["gene"], how="left")
        print(name)
    all_classes.columns = ['gene'] + samples.tolist()
    most_frequent_label_balanced = pd.DataFrame(all_classes.loc[:,'gene'])
    print('Pre-computing labels')
    for i, sample_out in enumerate(samples.tolist()):
        included = all_classes.loc[:, all_classes.columns != sample_out].mode(axis=1)[0]
        ## Get most frequent label:
        most_frequent_label_balanced[sample_out] = included
        print(i, sample_out)
    
    all_classes.to_csv(in_file_class_bal, index=False)
    most_frequent_label_balanced.to_csv(in_file_label_bal, index=False)
else:
    all_classes = pd.read_csv(in_file_class_bal)
    most_frequent_label_balanced = pd.read_csv(in_file_label_bal)

# Load UNBALANCED data
if not os.path.exists(in_file_class_unb) or not os.path.exists(in_file_label_unb):
    for i, name in enumerate(samples.tolist()):
        dev_m = pd.read_table(f'{path_in_unbalanced}/{name}_Dev_M.txt', names=colnames)[['gene', 'category']]
        tst_m = pd.read_table(f'{path_in_unbalanced}/{name}_Test_M.txt', names=colnames)[['gene', 'category']]
        trn_m = pd.read_table(f'{path_in_unbalanced}/{name}_Train_M.txt', names=colnames)[['gene', 'category']]
        all = pd.concat([dev_m, tst_m, trn_m])
        all = all[~all['gene'].isin(noisy_genes)]
        all_classes_unbalanced = all if i == 0 else pd.merge(all_classes_unbalanced, all, on=["gene"], how="left")
        print(name)
    
    all_classes_unbalanced.columns = ['gene'] + samples.tolist()
    most_frequent_label_unbalanced = pd.DataFrame(all_classes_unbalanced.loc[:,'gene'])
    print('Pre-computing labels')
    for i, sample_out in enumerate(samples.tolist()):
        included = all_classes_unbalanced.loc[:, all_classes_unbalanced.columns != sample_out].mode(axis=1)[0]
        ## Get most frequent label:
        most_frequent_label_unbalanced[sample_out] = included
        print(i, sample_out)
    
    all_classes_unbalanced.to_csv(in_file_class_unb, index=False)
    most_frequent_label_unbalanced.to_csv(in_file_label_unb, index=False)
else:
    all_classes_unbalanced = pd.read_csv(in_file_class_unb)
    most_frequent_label_unbalanced = pd.read_csv(in_file_label_unb)


# Calculate the acc
def calc_acc_per_sample(df_with_labels, samples):
    assert df_with_labels.shape[1] - 1 == len(samples)
    accuracies_per_sample = []
    for i, sample_out in enumerate(samples.tolist()):
        excluded = df_with_labels.loc[:, sample_out]
        included = df_with_labels.loc[:, df_with_labels.columns != sample_out]
        ## Get most frequent label:
        most_frequent_label = included.mode(axis=1)[0]
        current_sample_acc = pd.DataFrame({
            'most_frequent_label': most_frequent_label,
            'category': excluded.iloc[:,0]
        })
        ## Get accuracy in excluded sample
        current_sample_acc = current_sample_acc.query('most_frequent_label == category').shape[0] / current_sample_acc.shape[0]
        accuracies_per_sample.append(current_sample_acc)
        print(i, sample_out)
    return accuracies_per_sample

def calc_acc_per_sample_labels_precalc(df_with_labels, precalc_freq_labels, samples):
    assert df_with_labels.shape[1] - 1 == len(samples)
    accuracies_per_sample = []
    for i, sample_out in enumerate(samples.tolist()):
        ## Get most frequent label:
        current_sample_acc = pd.DataFrame({
            'most_frequent_label': precalc_freq_labels[sample_out],
            'category': df_with_labels[sample_out]
        })
        ## Get accuracy in excluded sample
        current_sample_acc = current_sample_acc.query('most_frequent_label == category').shape[0] / current_sample_acc.shape[0]
        accuracies_per_sample.append(current_sample_acc)
        print(i, sample_out)
    return accuracies_per_sample


accuracies_per_sample = calc_acc_per_sample_labels_precalc(all_classes, most_frequent_label_balanced, samples)
accuracies_per_sample_unbalanced = calc_acc_per_sample_labels_precalc(all_classes_unbalanced, most_frequent_label_unbalanced, samples)

plt.figure(figsize=(6, 10))
sns.boxplot(y=accuracies_per_sample, color='lightblue')
sns.swarmplot(y=accuracies_per_sample, color='darkblue')
plt.title('[Baseline Balanced Split]\nAccuracy of Label Prediction per Sample (N=49)')
plt.xlabel('Samples')
plt.ylabel('Accuracy')
plt.ylim(0.7,0.95)
plt.savefig('FILTERED_Baseline_category_accuracy_of_label_prediction_per_sample_balanced_split.pdf')
plt.close()

plt.figure(figsize=(6, 10))
sns.boxplot(y=accuracies_per_sample_unbalanced, color='lightblue')
sns.swarmplot(y=accuracies_per_sample_unbalanced, color='darkblue')
plt.title('[Baseline UNbalanced Split]\nAccuracy of Label Prediction per Sample (N=49)')
plt.xlabel('Samples')
plt.ylabel('Accuracy')
plt.ylim(0.7,0.95)
plt.savefig('FILTERED_Baseline_category_accuracy_of_label_prediction_per_sample_UNbalanced_split.pdf')
plt.close()


df = pd.DataFrame({'Balanced': accuracies_per_sample, 'Unbalanced': accuracies_per_sample_unbalanced})
t_stat, p_value = ttest_rel(accuracies_per_sample, accuracies_per_sample_unbalanced)
# Plot
sns.boxplot(data=df)
plt.title(f'Paired t-test\np-value: {p_value:.1e}')
# Add significance if p < 0.05
if p_value < 0.05:
    y, h, col = df.max().max() + 0.01, 0.01, 'k'
    plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text(0.5, y+h, "*", ha='center', va='bottom', color=col)

plt.show()
plt.savefig('FILTERED_Baseline_category_accuracy_of_label_prediction_per_sample_paired.pdf')


# __IMPORTANT__
# For this section below run first 
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_variable_genes/fcdnn_variable_genes.sh
# __^^^^^^^^^__

variable_genes = pd.read_csv('/tmp/all_any_variable_genes.txt', names = ['gene'])

all_classes.reset_index(inplace=True)
all_classes_unbalanced.reset_index(inplace=True)

most_frequent_label_balanced.reset_index(inplace=True)
most_frequent_label_unbalanced.reset_index(inplace=True)

# Get the mayority vote labes of most variable genes:
variable_genes_indexes = all_classes.reset_index().merge(variable_genes,how='inner', on=['gene'])['index']
all_classes = all_classes[all_classes['index'].isin(variable_genes_indexes)].drop(columns=['index'])
all_classes_unbalanced = all_classes_unbalanced[all_classes_unbalanced['index'].isin(variable_genes_indexes)].drop(columns=['index'])
most_frequent_label_balanced = most_frequent_label_balanced[most_frequent_label_balanced['index'].isin(variable_genes_indexes)].drop(columns=['index'])
most_frequent_label_unbalanced = most_frequent_label_unbalanced[most_frequent_label_unbalanced['index'].isin(variable_genes_indexes)].drop(columns=['index'])

var_genes_accuracies_per_sample = calc_acc_per_sample_labels_precalc(all_classes, most_frequent_label_balanced, samples)
var_genes_accuracies_per_sample_unbalanced = calc_acc_per_sample_labels_precalc(all_classes_unbalanced, most_frequent_label_unbalanced, samples)

# Plot
df = pd.DataFrame({'Balanced': var_genes_accuracies_per_sample, 'Unbalanced': var_genes_accuracies_per_sample_unbalanced})
t_stat, p_value = ttest_rel(var_genes_accuracies_per_sample, var_genes_accuracies_per_sample_unbalanced)
# Plot
sns.boxplot(data=df)
sns.stripplot(data=df, color='black')
plt.title(f'Paired t-test\np-value: {p_value:.1e}')
# Add significance if p < 0.05
if p_value < 0.05:
    y, h, col = df.max().max() + 0.01, 0.01, 'k'
    plt.plot([0, 0, 1, 1], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text(0.5, y+h, "*", ha='center', va='bottom', color=col)

plt.show()
plt.savefig('FILTERED_Baseline_category_accuracy_of_label_prediction_per_sample_paired_variable_genes.pdf')
plt.close()

