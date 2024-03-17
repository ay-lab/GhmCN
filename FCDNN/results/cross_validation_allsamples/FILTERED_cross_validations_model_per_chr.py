# srun --nodes=1 --ntasks=1 --cpus-per-task=24 --mem=80g --time=04:00:00 --pty bash -i
# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/code/; python

# NOTE:
# 1c part 2: The reviewers specifically asked us to do some robustness analysis where we train and test
# multiple models on different combinations of chromosomes. You did this for 10 random chrs (holding one 
# chr off the training at a time) as far as I remember. My suggestion was to do this for each chr simply 
# so you train 19 different models where for each model you train on 17 chrs, tune on 1 and test on 1 
# (like train 1:17, tune 18, test 19) etc. And then the per chr plot was going to be made with those 
# different models where for each chr you would have 49 dots (accuracy or AUC of predictions for that 
# chromosome for each sample from the model that specific chr was held out). So 19 different models, 
# 49 dots per each chromosome corresponding to 49 samples. Let me know if this is not clear. We have to clarify.

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
    return x_train, x_test, x_dev, y_train_hot, y_test, y_test_hot, y_dev_hot, genes_in_chrom


# If first time, then True, if not ,False and jump to analysis:
compress_data=False

if compress_data:
    # Load and save as binary numpy (X,Y) Train, Test, Dev
    x_data, y_data, y_data_hot = load_dataset("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples", "Test", 1)
    x_data_dev, y_data_dev, y_data_hot_dev = load_dataset("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples", "Dev", 1)
    x_data_train, y_data_train, y_data_hot_train = load_dataset("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples", "Train", 1)
    
    x_all_data = np.hstack((np.hstack((x_data, x_data_dev)), x_data_train))
    y_all_data = np.concatenate((np.concatenate((y_data, y_data_dev)), y_data_train))
    
    x_all_data.shape
    y_all_data.shape
    np.save('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_x.npy', x_all_data)
    np.save('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_y.npy', y_all_data)
    
    # Load and save as binary pandas (M) Train, Test, Dev
    m_data = load_m_dataset("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples", "Test")
    m_data_dev = load_m_dataset("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples", "Dev")
    m_data_train = load_m_dataset("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples", "Train")
    m_data = pd.concat([m_data, m_data_dev, m_data_train]).reset_index(drop=True)
    m_data.to_parquet('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_m.parquet')
    del x_data, y_data, y_data_hot, x_data_dev, y_data_dev, y_data_hot_dev, x_data_train, y_data_train, y_data_hot_train, x_all_data, y_all_data, m_data, m_data_dev, m_data_train

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
samples_in_dataset = m_all_data.groupby('gene').size().mean().astype(int)

# Pre-calculate genes per chrom
chrom_genes_count = m_all_data.groupby('chr')['gene'].count()//samples_in_dataset

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

# Save the dict indicating the chr excluded for test and dev
with open(output_dict_file, 'w') as f:
    json.dump(chr_test_dev, f)

results=[]
# Run throught a single chromomosome (sys.args) or the whole series
chromosomes = [unique_chromosomes[subset_chrom]] if subset_chrom is not None else unique_chromosomes
for chrom in chromosomes:
    # Declare outputs
    print(f'\n# >>> training {chrom} model')
    
    output_figure_file = os.path.join(output_dir,f'FILTERED_accuracy_plot_cross_validations_{chrom}.png')
    output_model_file = os.path.join(output_dir,f'fcdnn_model_parameters_{chrom}.pkl')
    output_chrom = os.path.join(output_dir,f'FILTERED_cross_validation_all_scores_per_sample_{chrom}.tsv')
    
    x_train, x_test, x_dev, y_train_hot, y_test, y_test_hot, y_dev_hot, genes_in_chrom = subset_dataset(m_all_data, chr_masks, chrom_genes_count, chrom)
    parameters, _, _, _ = FCDNN(x_train, y_train_hot, x_dev, y_dev_hot, layers_dim = layers,
                                Parameters=None, learning_rate = 0.0001,
                                num_epochs = 60,  minibatch_size = 512, generate_plots = True,
                                AppendEach = 5, keepProb=0.85, betaReg=0.01, seed = 1921,print_acc = False,
                                decay = 0.975, decay_schedule = 30, OutputFigure=output_figure_file)
    
    # Save model parameters
    with open(output_model_file, 'wb') as f:
        pickle.dump(parameters, f, protocol=4)
    
    # calculate chromosome-wide metrics
    obs, pred, prob, acc, _ = run_predictions_no_save(x_test, y_test, y_test_hot, parameters, output_layer)
    _, _, auc_score, f1 = calculate_metrics(obs, pred, prob, acc, output_layer, print_metrics=False)
    results.append((chrom,'all_samples',acc,round(auc_score,4), round(f1,4)))
    
    # calculate unbiased metrics
    print(f'\n# >>> calculating {chrom} unbiased metrics per sample\n')
    for i in range(samples_in_dataset):  # +1 to include the last chunk that may be smaller than 1340
        sample_name = samples[i]
        start_index = i * genes_in_chrom
        end_index = start_index + genes_in_chrom
        # Slicing the arrays to get the current chunk
        obs_subset = obs[start_index:end_index]
        pred_subset = pred[start_index:end_index]
        prob_subset = prob[start_index:end_index]
        
        # Ensure we're not processing an empty chunk at the end
        if len(obs_subset) > 0:
            acc = round(sum(obs_subset == pred_subset)/len(obs_subset),4)
            # Calculate metrics for the current chunk
            _, _, auc_score, f1 = calculate_metrics(obs_subset, pred_subset, prob_subset, acc, output_layer, print_metrics=False)
            # Append the results to the list
            results.append((chrom,sample_name,acc, round(auc_score,4), round(f1,4)))
    
    # Save per chrom
    merged_results = pd.DataFrame(results, columns = ['chr','sample','accuracy','auroc','f1'])
    merged_results.to_csv(output_chrom,index=False)            

# Save scores
output_all = os.path.join(output_dir,'FILTERED_cross_validation_all_scores_per_sample_per_chrom.tsv')
merged_results = pd.DataFrame(results, columns = ['chr','sample','accuracy','auroc','f1'])
merged_results.to_csv(output_all,index=False)
