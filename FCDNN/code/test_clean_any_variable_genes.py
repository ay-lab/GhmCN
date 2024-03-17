# Sconda conda activate TensorFlow_CPU_01; cd /mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/FCDNN/code; python
#!/mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/envs/TensorFlow_CPU_01/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import tensorflow as tf
from pathlib import Path
from tf_utils import *

from typing import List
from sklearn.metrics import precision_score, recall_score, f1_score, roc_curve, auc

data_path = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/'
out_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_variable_genes/'
model_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/rerun/C2B/training/merged/models/'
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/rerun/C2B/training/merged/models/fcdnn_model_parameters_CleanSamples.pkl
name="CleanSamplesAnyVariable"
data_prefix = os.path.join(data_path, name)
categories=1
suffix_x = "_X.npz"
suffix_y = "_Y.npz"
x_data = np.load(data_prefix + suffix_x)['arr_0']
y_data = np.load(data_prefix + suffix_y)['arr_0']
y_data_hot = one_hot_matrix(y_data, categories) if categories > 1 else y_data
print ('X shape:\t%s\nY shape:\t%s' % ( x_data.shape, y_data_hot.shape))

# calculate accuracy
def run_predictions(x_data, y_data, y_data_hot, parameters, output_layer, output_file_pred):
    prob, pred, acc = predict(x_data, y_data_hot, parameters, output_layer)
    highest_prob = np.max(prob,axis = 0)
    pred = pred.flatten() if output_layer == 'sigmoid' else pred
    obs = np.squeeze(y_data.astype(int))
    side2prob = pd.DataFrame({'Observed': obs,'Predicted': pred, 'Probability': highest_prob})
    side2prob.to_csv(output_file_pred, index=False)
    return obs, pred, highest_prob, acc, prob

# calculate unbiased metrics
def calculate_metrics(obs, pred, prob, acc, output_layer, output_file_metrics = None, name = None, save_metrics:bool=False):
    if output_layer == "sigmoid":
        precision = precision_score(obs, pred)
        recall = recall_score(obs, pred)
        f1 = f1_score(obs, pred)
        fpr, tpr, _ = roc_curve(obs, prob)
        auc_score = auc(fpr, tpr)
        print("Precision:\t%.4f\nRecall: \t%.4f\nF1 score:\t%.4f\nAUC score:\t%.4f"  %  (precision, recall, f1, auc_score) )
    else:
        precision = precision_score(obs, pred, average = "weighted")
        recall = recall_score(obs, pred, average = "weighted")
        f1 = f1_score(obs, pred, average = "weighted")
        fpr, tpr, _ = np.nan, np.nan, np.nan
        auc_score = np.nan
        print("Precision:\t%.4f\nRecall: \t%.4f\nF1 score:\t%.4f\nAUC score:\t%.4f"  %  (precision, recall, f1, auc_score) )
    if save_metrics:
        metrics = pd.DataFrame({'Accuracy': acc,
                                'Precision' : precision,
                                'Recall' : recall,
                                'F1_score' : f1,
                                'AUC_score' : auc_score}, index=[name]).round(4)
        metrics.to_csv(output_file_metrics)
    return fpr, tpr, auc_score, f1

# ptotting
def plot_roc(fpr, tpr, auc_score, f1, name, output_file_roc):
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.4f)' % auc_score)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f"Receiver Operating Characteristic (ROC) Curve\n{name}")
    plt.legend(loc='lower right')
    plt.text(0.825, 0.125, f'F1: {f1:.4f}', transform=plt.gca().transAxes)
    plt.savefig(output_file_roc)
    plt.close()


metrics_file_name = 'metrics_CleanSamples_tested_on_AnyVariableGenes.csv'
predictions_file_name = 'predictions_CleanSamples_tested_on_AnyVariableGenes.csv'
output_file_roc = os.path.join(out_path,'roc_curve_CleanSamples_tested_on_AnyVariableGenes.png')
# declare parameters to load
parameters_file_model = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/rerun/C2B/training/merged/models/fcdnn_model_parameters_CleanSamples.pkl'
# declate metrics outputs
output_file_metrics = os.path.join(out_path,metrics_file_name)
# declate predictions output
output_file_pred = os.path.join(out_path,predictions_file_name)

# load parameters
with open(parameters_file_model, 'rb') as f:
    parameters = pickle.load(f)

# calculate accuracy
obs, pred, prob, acc, _ = run_predictions(x_data, y_data, y_data_hot, parameters, 'sigmoid', output_file_pred)
output_file_pred
# calculate unbiased metrics
fpr, tpr, auc_score, f1 = calculate_metrics(obs, pred, prob, acc, 'sigmoid', output_file_metrics, name, True)
# Accuracy:       0.7838
# Precision:      0.7849
# Recall:         0.7818
# F1 score:       0.7833
# AUC score:      0.8627