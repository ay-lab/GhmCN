# srun --nodes=1 --ntasks=1 --cpus-per-task=24 --mem=80g --time=04:00:00 --pty bash -i
# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/code/; python

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import re
import json

import sys
sys.path.append('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/code')
from tf_utils import *
from fcdnn import *

# Convert the first argument to an integer
try:
    subset_chrom = int(sys.argv[1])
except:
    subset_chrom = None


def load_m_dataset(data_prefix: str, subset: str):
    if subset not in ["Train", "Dev", "Test"]:
        raise ValueError("Invalid value for 'subset'. Allowed values are 'Train', 'Dev', and 'Test'.")
    suffix_m = "_"+subset+"_M.txt"
    m_data = pd.read_table(data_prefix + suffix_m, names = ["gene","category","binary_category","chr", "order","tpm"])
    return m_data


def run_predictions_no_save(x_data, y_data, y_data_hot, parameters, output_layer):
    prob, pred, acc = predict(x_data, y_data_hot, parameters, output_layer)
    highest_prob = np.max(prob,axis = 0)
    pred = pred.flatten() if output_layer == 'sigmoid' else pred
    obs = y_data.astype(int)
    return obs, pred, highest_prob, acc, prob


def subset_dataset(m_all_data, chr_masks, chrom_genes_count, chrom):
    # Use pre-calculated boolean masks
    test_idx = m_all_data[chr_masks['testing'][chrom]].index
    dev_idx = m_all_data[chr_masks['dev'][chrom]].index
    train_idx = m_all_data[chr_masks['train'][chrom]].index
    genes_in_test = m_all_data[chr_masks['testing'][chrom]].reset_index(drop=True)
    genes_in_chrom = chrom_genes_count[chrom]
    
    # Subset training
    x_train =  x_all_data[:, train_idx]
    y_train =  y_all_data[train_idx]
    # Subset test
    x_test =  x_all_data[:, test_idx]
    y_test =  y_all_data[test_idx]
    # Subset dev
    x_dev =  x_all_data[:, dev_idx]
    y_dev =  y_all_data[dev_idx]
    # reshape
    y_train_hot = y_train.reshape((1, y_train.shape[0]))
    y_test_hot = y_test.reshape((1, y_test.shape[0]))
    y_dev_hot = y_dev.reshape((1, y_dev.shape[0]))
    return x_train, x_test, x_dev, y_train_hot, y_test, y_test_hot, y_dev_hot, genes_in_chrom, genes_in_test


# Load data:
x_all_data = np.load('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_x.npy')
y_all_data = np.load('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_y.npy')
m_all_data = pd.read_parquet('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_m.parquet')
samples = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/filtered_files_surnames.csv', header=None)[7].tolist()

# constants
output_dir = "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/model_per_chr/"
output_dict_file = os.path.join(output_dir,'chr_excluded_for_test_and_dev.json')
cross_validations_auc_rocs = []
cross_validations_f1s = []
chrom_order = []
# categories, layers
categories = 1 # Total labels: set to 1 if binary.
layers = [230,200,100,50,1]


# Subest 10 chromosomes:
unique_chromosomes = m_all_data['chr'].unique()
unique_chromosomes = unique_chromosomes[(unique_chromosomes != "chrY") & (unique_chromosomes != "chrX")]
output_layer = "sigmoid"

# Pre-calculate number of samples composing dataset
n_samples = m_all_data.groupby('gene').size().mean().astype(int)

# Pre-calculate genes per chrom
chrom_genes_count = m_all_data.groupby('chr')['gene'].count()//n_samples

# Pre-calculate boolean masks for each chromosome to avoid repetitive indexing
chr_masks = {'testing':{}, 'dev':{}, 'train':{}}
chr_test_dev = {}
for i, chrom in enumerate(unique_chromosomes):
    # Find the next chromosome in the list, wrapping around to the first if at the end
    next_chrom = unique_chromosomes[(i + 1) % len(unique_chromosomes)]
    chr_test_dev[chrom] = next_chrom
    # build dict of dict of masks per test chrom
    chr_masks['testing'][chrom] = m_all_data["chr"] == chrom
    chr_masks['dev'][chrom] = m_all_data["chr"] == next_chrom
    chr_masks['train'][chrom] = ~(
        (m_all_data["chr"] == chrom) | (m_all_data["chr"] == next_chrom)
    )

# Load the variable genes file
variable_genes = pd.read_csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences/FILTERED_any_degree_of_variable_label_genes_promoter_coords.bed", sep="\t", header=None)[[0,3]]
variable_genes.columns = ['chr','gene']

# initialize dict of results per sample:
results = {}
for i in range(n_samples):
    results[i] = {'obs': [], 'pred': []}

# Save the dict indicating the chr excluded for test and dev
# Run throught a single chromomosome (sys.args) or the whole series
output_dir = "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/model_per_chr/"
for chrom in unique_chromosomes:
    print(chrom)
    output_model_file = os.path.join(output_dir,f'fcdnn_model_parameters_{chrom}.pkl')
    # Get model from chr
    with open(output_model_file, 'rb') as f:
        parameters = pickle.load(f)
    
    # Get data from chr
    x_train, x_test, x_dev, y_train_hot, y_test, y_test_hot, y_dev_hot, genes_in_chrom, genes_in_test = subset_dataset(m_all_data, chr_masks, chrom_genes_count, chrom)
    # Get variable genes from chr
    chr_var_genes = variable_genes.query('chr == @chrom').gene.tolist()
    # Get position of variable genes
    subset_predictions = genes_in_test[genes_in_test.gene.isin(chr_var_genes)].index.tolist()
    n_genes = len(subset_predictions)//n_samples
    # obtain observed and predicted
    x_test_var_genes = x_test[:,subset_predictions]
    y_test_hot_var_genes = y_test_hot[:,subset_predictions]
    _, pred, _ = predict(x_test_var_genes, y_test_hot_var_genes, parameters, output_layer)
    obs = y_test[subset_predictions]
    pred = pred.flatten()
    for i in range(n_samples):
        idx_start = n_genes*i
        idx_end = idx_start+n_genes
        # Append obs and pred directly to the lists
        results[i]['obs'].append(obs[idx_start:idx_end])
        results[i]['pred'].append(pred[idx_start:idx_end])

# Calculate real acc in variable samples:
acc_per_sample = []
for i in range(n_samples):
    obs = np.concatenate(results[i]['obs'])
    pred = np.concatenate(results[i]['pred'])
    acc = np.sum(obs == pred) / len(obs)
    acc_per_sample.append(acc)

pd.DataFrame({"sample_number": range(n_samples), "acc":acc_per_sample}).to_csv(output_dir+'/corrected_accuracies_in_variable_genes.csv', index=False)