# `CleanSamples` model
# srun --nodes=1 --ntasks=1 --cpus-per-task=24 --mem=100g --time=04:00:00 --pty bash -i
# Sconda conda activate  tensorflow2; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples; python

import json
import re
import matplotlib.pyplot as plt
import sys
import os
import pkg_resources
import numpy as np
import pandas as pd
import seaborn as sns

from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    roc_curve,
    roc_auc_score,
)
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import (
    Dense,
    Dropout,
)
from tensorflow.keras.initializers import GlorotUniform

from tensorflow.keras.optimizers import Adam
from keras.models import model_from_json


## Declare some constants
cross_validations_auc_rocs = []
cross_validations_f1s = []
chrom_order = []
# categories, layers
categories = 1  # Total labels: set to 1 if binary.
layers = [230, 200, 100, 50, 1]

# output figures:
fig_out_dir = "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/keras_models_figures"
model_out_dir = "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/keras_models/"
file_path = os.path.join(
    fig_out_dir, "FILTERED_cross_validation_all_scores_per_sample_per_chrom.tsv"
)

### Loading Datasets
x_all_data = np.load(
    "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_x.npy"
).T
y_all_data = np.load(
    "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_y.npy"
).reshape(-1, 1)
m_all_data = pd.read_parquet(
    "/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_m.parquet"
)

print(
    f"""
X all data shape:\t{x_all_data.shape}
Y all data shape:\t{y_all_data.shape}
Meta data shape:\t{m_all_data.shape}
"""
)
### Funciton to remove layer that is causing issues:
# Load the original JSON data
def json_remove_input_layer(json_file):
    # Read json
    with open(json_file, "r") as file:
        data = json.load(file)
        file.close()
    # Check if it has the problem and delete
    if (
        data["config"]["layers"]
        and data["config"]["layers"][0]["class_name"] == "InputLayer"
    ):
        del data["config"]["layers"][0]
        print("deleted 'InputLayer' layer config")
    
    # Save the modified data back to a new JSON file
    with open(json_file, "w") as file:
        json.dump(data, file, indent=4)
        file.close()


# Remove the first layer from the 'layers' list within the 'config' key
# Subest 10 chromosomes:
unique_chromosomes = m_all_data["chr"].unique()
unique_chromosomes = unique_chromosomes[
    (unique_chromosomes != "chrY") & (unique_chromosomes != "chrX")
]
output_layer = "sigmoid"


# Pre-calculate boolean masks for each chromosome to avoid repetitive indexing
chr_masks_testing = {chrom: m_all_data["chr"] == chrom for chrom in unique_chromosomes}
# Create the chr_masks_training dictionary
chr_masks_testing = {}
chr_masks_dev = {}
chr_masks_training = {}
for i, chrom in enumerate(unique_chromosomes):
    next_chrom = unique_chromosomes[(i + 1) % len(unique_chromosomes)]
    chr_masks_testing[chrom] = m_all_data["chr"] == chrom
    chr_masks_dev[chrom] = m_all_data["chr"] == next_chrom
    # Find the next chromosome in the list, wrapping around to the first if at the end
    chr_masks_training[chrom] = ~(
        (m_all_data["chr"] == chrom) | (m_all_data["chr"] == next_chrom)
    )


# Pre-calculate number of samples composing dataset
samples_in_dataset = m_all_data.groupby("gene").size().mean().astype(int)

# Pre-calculate genes per chrom
chrom_genes_count = m_all_data.groupby("chr")["gene"].count() // samples_in_dataset
pd.DataFrame(chrom_genes_count).reset_index().to_csv(f'{model_out_dir}/genes_per_chrom.csv', index=False)
# Make this a function to get a cleaner code:
def make_keras_model(X):
    keras_model = Sequential()
    keras_model.add(
        Dense(
            200,
            activation="relu",
            name="1",
            input_dim=X.shape[1],
            kernel_initializer=GlorotUniform(),
        )
    )
    keras_model.add(Dropout(0.15))
    keras_model.add(
        Dense(100, activation="relu", name="2", kernel_initializer=GlorotUniform())
    )
    keras_model.add(Dropout(0.15))
    keras_model.add(
        Dense(50, activation="relu", name="3", kernel_initializer=GlorotUniform())
    )
    keras_model.add(
        Dense(
            1, activation="sigmoid", name="Output", kernel_initializer=GlorotUniform()
        )
    )
    keras_model.compile(
        optimizer=Adam(learning_rate=0.0001),
        loss="binary_crossentropy",
        metrics=["accuracy"],
    )
    return keras_model


def test_and_auc(keras_model, X, Y, samples_in_dataset, genes_in_chrom, chrom):
    results = []
    for i in range(
        samples_in_dataset
    ):  # +1 to include the last chunk that may be smaller than 1340
        start_index = i * genes_in_chrom
        end_index = start_index + genes_in_chrom
        # Slicing the arrays to get the current chunk
        X_sub = X[start_index:end_index, :]
        Y_sub = Y[start_index:end_index]
        preds = keras_model.evaluate(x=X_sub, y=Y_sub)
        ProbableValues = keras_model.predict(X_sub)
        Obs = pd.Series(Y_sub.flatten().astype(int))
        Pred = pd.Series(ProbableValues.round().flatten().astype(int))
        auc_roc = roc_auc_score(Obs, Pred)
        results.append((chrom, f"sample_{i}", round(auc_roc, 4), round(preds[1], 4)))
    
    return results


### Train Model
results = []
for chrom in unique_chromosomes:
    # Declare outputs
    output_model_configs = os.path.join(model_out_dir, "keras_model_" + chrom + ".json")
    output_model_weights = os.path.join(model_out_dir, "keras_model_" + chrom + ".h5")

    # Use pre-calculated boolean masks
    test_idx = m_all_data[chr_masks_testing[chrom]].index
    dev_idx = m_all_data[chr_masks_dev[chrom]].index
    train_idx = m_all_data[chr_masks_training[chrom]].index
    genes_in_chrom = chrom_genes_count[chrom]

    # Subset the x and y datasets using the obtained indices
    X = x_all_data[train_idx, :]
    Y = y_all_data[train_idx]
    x_dev = x_all_data[dev_idx, :]
    y_dev = y_all_data[dev_idx]
    x_test = x_all_data[test_idx, :]
    y_test = y_all_data[test_idx]

    # define the keras model
    clean_data_model = make_keras_model(X)
    clean_data_model.fit(x=X, y=Y, batch_size=1024, epochs=60)
    preds = clean_data_model.evaluate(x=x_test, y=y_test)

    # serialize model to JSON
    model_json = clean_data_model.to_json()
    with open(output_model_configs, "w") as json_file:
        json_file.write(model_json)
        json_file.close()
    # Load back for formatting
    json_remove_input_layer(output_model_configs)
    # serialize weights to HDF5
    clean_data_model.save_weights(output_model_weights)

    results.append(
        test_and_auc(
            clean_data_model, x_test, y_test, samples_in_dataset, genes_in_chrom, chrom
        )
    )


flat_results = [item for sublist in results for item in sublist]
flat_results = pd.DataFrame(
    flat_results, columns=["chr_tested_on", "sample", "auroc", "accuracy"]
)
flat_results.to_csv(file_path, index=False)
# Make plots
# Load the CSV file into a DataFrame
auc_data = pd.read_csv(file_path)


mean_auc_score = 0.869


# Re-defining the sort function with the 're' module imported
def sort_chromosomes(chromosomes):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(chromosomes, key=alphanum_key)


# Applying the custom sort function to the chromosome names again
sorted_chromosomes = sort_chromosomes(auc_data["chr_tested_on"].unique())

# Now, creating the boxplot chromosome-wise
plt.figure(figsize=(16, 10))
sns.set(style="whitegrid")

# Creating the boxplot with chromosomes sorted correctly
sns.stripplot(
    x="chr_tested_on",
    y="auroc",
    data=auc_data,
    color="black",
    order=sorted_chromosomes,
    size=4,
    jitter=True,
    alpha=0.6,
)
ax = sns.boxplot(
    x="chr_tested_on",
    y="auroc",
    data=auc_data,
    color="lightblue",
    order=sorted_chromosomes,
)

# Adding black points for each data point on top of the boxplots

# Adding a red line for the mean AUC score
plt.axhline(
    mean_auc_score,
    color="red",
    linestyle="--",
    linewidth=1,
    label=f"Combined Model mean AUROC: {mean_auc_score:.2f}",
)

# Adding a legend
plt.legend(loc="upper right", fontsize=14, frameon=True)

# Customizing the plot
ax.set_title("AUC Distribution per Chromosomes", fontsize=20)
ax.set_xlabel("Chromosome", fontsize=16)
ax.set_ylabel("AUC Value", fontsize=16)
ax.tick_params(axis="x", labelsize=14)
ax.tick_params(axis="y", labelsize=14)

plt.tight_layout()

# Save the figure with the mean AUC line
output_sample_boxplot_path = os.path.join(fig_out_dir,
    "FILTERED_cross_validation_auc_distribution_per_chrom_boxplot_overlayed"
)
plt.savefig(output_sample_boxplot_path + ".png", dpi=300, bbox_inches="tight")
plt.savefig(output_sample_boxplot_path + ".pdf", dpi=150, bbox_inches="tight")

# Show the plot
plt.show()
plt.close()