# srun --nodes=1 --ntasks=1 --cpus-per-task=24 --mem=80g --time=04:00:00 --pty bash -i
# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/code/; python
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import csv
import random

import sys
sys.path.append('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/code')
from tf_utils import *
from fcdnn import *

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

cross_validations_auc_rocs = []
cross_validations_f1s = []
chrom_order = []
# categories, layers
categories = 1 # Total labels: set to 1 if binary.
layers = [230,200,100,50,1]

# output figures:
output_dir = "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples"

# Subest 10 chromosomes:
unique_chromosomes = m_all_data['chr'].unique()
unique_chromosomes = unique_chromosomes[(unique_chromosomes != "chrY") & (unique_chromosomes != "chrX")]
output_layer = "sigmoid"

# Pre-calculate boolean masks for each chromosome to avoid repetitive indexing
chr_masks = {chrom: m_all_data['chr'] == chrom for chrom in unique_chromosomes}

# Pre-calculate number of samples composing dataset
samples_in_dataset = m_all_data.groupby('gene').size().mean().astype(int)

# Pre-calculate genes per chrom
chrom_genes_count = m_all_data.groupby('chr')['gene'].count()//samples_in_dataset

results=[]
# Prepare dictionaries to hold the subsets
for chrom in unique_chromosomes:
    # Declare outputs
    output_figure_file = os.path.join(output_dir,'FILTERED_accuracy_plot_cross_validations_'+chrom+'.png')
    
    # Use pre-calculated boolean masks
    test_idx = m_all_data[chr_masks[chrom]].index
    train_idx = m_all_data[~chr_masks[chrom]].index
    genes_in_chrom = chrom_genes_count[chrom]
    
    # Subset the x and y datasets using the obtained indices
    x_train =  x_all_data[:, train_idx]
    y_train =  y_all_data[train_idx]
    x_test =  x_all_data[:, test_idx]
    y_test =  y_all_data[test_idx]
    
    y_train_hot = y_train.reshape((1, y_train.shape[0]))
    y_test_hot = y_test.reshape((1, y_test.shape[0]))
    
    parameters, _, _, _ = FCDNN(x_train, y_train_hot, x_test, y_test_hot, layers_dim = layers,
                                Parameters=None, learning_rate = 0.0001,
                                num_epochs = 6,  minibatch_size = 64, generate_plots = True,
                                AppendEach = 5, keepProb=0.85, betaReg=0.01, seed = 1921,
                                decay = 0.975, decay_schedule = 30, OutputFigure=output_figure_file)
    
    obs, pred, prob, acc, _ = run_predictions_no_save(x_test, y_test, y_test_hot, parameters, output_layer)
    
    # calculate unbiased metrics
    for i in range(samples_in_dataset + 1):  # +1 to include the last chunk that may be smaller than 1340
        start_index = i * genes_in_chrom
        end_index = start_index + genes_in_chrom
        # Slicing the arrays to get the current chunk
        obs_subset = obs[start_index:end_index]
        pred_subset = pred[start_index:end_index]
        prob_subset = prob[start_index:end_index]
        
        # Ensure we're not processing an empty chunk at the end
        if len(obs_subset) > 0:
            # Calculate metrics for the current chunk
            _, _, auc_score, f1 = calculate_metrics(obs_subset, pred_subset, prob_subset, acc, output_layer, print_metrics=False)
            # Append the results to the list
            results.append((chrom,f'sample_{i}',round(auc_score,4), round(f1,4)))
    
    _, _, auc_score, f1 = calculate_metrics(obs, pred, prob, acc, output_layer)
    
    cross_validations_auc_rocs.append(round(auc_score,4))
    cross_validations_f1s.append(round(f1,4))
    chrom_order.append(chrom)

# Save scores
output_chrom = os.path.join(output_dir,'FILTERED_chrom_order.tsv')
output_auc = os.path.join(output_dir,'FILTERED_auc_scores_cross_validation.tsv')
output_f1 = os.path.join(output_dir,'FILTERED_f1_scores_cross_validation.tsv')
output_all = os.path.join(output_dir,'FILTERED_cross_validation_all_scores_per_sample_per_chrom.tsv')
results = pd.DataFrame(results, columns = ['chr','sample','auc','f1'])
results.to_csv(output_all,index=False)

with open(output_chrom, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')  
    for score in chrom_order:
        writer.writerow([score])

with open(output_auc, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')  
    for score in cross_validations_auc_rocs:
        writer.writerow([score])

with open(output_f1, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')  
    for score in cross_validations_f1s:
        writer.writerow([score])

# Make plot of scores and save scores
plt.boxplot(cross_validations_auc_rocs)
output_fig = os.path.join(output_dir,'FILTERED_auc_scores_cross_validation.png')
plt.title('Boxplot of AUC Cross Validations')
plt.ylabel('AUC Scores')
plt.savefig(output_fig)
plt.close()

plt.boxplot(cross_validations_f1s)
output_fig = os.path.join(output_dir,'FILTERED_f1_scores_cross_validation.png')
plt.title('Boxplot of F1 from Cross Validations')
plt.ylabel('F1 Scores')
plt.savefig(output_fig)
plt.close()
