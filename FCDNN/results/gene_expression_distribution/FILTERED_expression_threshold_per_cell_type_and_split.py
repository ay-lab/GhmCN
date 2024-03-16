# Sconda conda activate abc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/gene_expression_class_concordance; python
# python expression_threshold_per_celltype.py
import pandas as pd
import numpy as np
import os


# Function to TPM normalize the gene counts
def tpm_normalize(df):
    # Assuming column 6 has gene length in base pairs (bp)
    length_kb = df.iloc[:, 5] / 1000  # Convert gene length from bp to kilobases (kb)
    # Assuming column 7 has raw read counts
    reads = df.iloc[:, 6]
    # Calculate Reads Per Kilobase (RPK)
    rpk = reads / length_kb
    # Calculate the scaling factor by summing all RPK values and dividing by 1,000,000
    scaling_factor = rpk.sum() / 1_000_000
    # Divide RPK by the scaling factor to get TPM
    tpm = rpk / scaling_factor
    return tpm


# Get the gene expression thresholds for balanced and biological splits:
def get_expression_thresholds(df):
    tpm_values = tpm_normalize(df)
    non_zero_genes = tpm_values[tpm_values > 0]
    bio_threshold = np.median(non_zero_genes)
    balanced_threshold = np.median(tpm_values)
    return bio_threshold, balanced_threshold


# Read the CSV file
base_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/gene_expression_class_concordance'
out_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/gene_expression_distribution'
out_file=os.path.join(out_path,'FILTERED_median_threshold_per_cell_type.csv')

samples = pd.read_csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/filtered_files_surnames.csv", header = None, usecols = [7,10])
samples.columns = ["sample_name","rna_file"]

for index, row in samples.iterrows():
    bio_thresholds = []
    balanced_thresholds = []
    # Initialize lists to store labels and concordance values
    # Load and process each replicate
    rna_data = pd.read_csv(row['rna_file'], sep='\t', comment='#')
    bio, balanced = get_expression_thresholds(rna_data)
    bio_thresholds.append(bio)
    balanced_thresholds.append(balanced)
    # Save median thresholds per cell type
    result_df = pd.DataFrame({'cell_type': row['sample_name'],
                            'bio_threshold_tpm': [np.mean(bio_thresholds)],
                            'balanced_threshold_tpm': [np.mean(balanced_thresholds)]})
    result_df.to_csv(out_file, mode='a', header=(not os.path.exists(out_file)), index=False)

# Calculate the median and +/- 95% values
data = pd.read_csv(out_file)
median_bio_threshold = data.iloc[:, 1].median()
median_balanced_threshold = data.iloc[:, 2].median()

# Calculate the 95% confidence interval for the median
# Since we can't calculate a true confidence interval for a median without making assumptions
# about the underlying distribution, we'll use the percentile approach as an approximation.

# Calculate the 2.5th and 97.5th percentiles for the second and third columns
percentile_2_5_bio_threshold = data.iloc[:, 1].quantile(0.025)
percentile_2_5_balanced_threshold = data.iloc[:, 2].quantile(0.025)
percentile_25_bio_threshold = data.iloc[:, 1].quantile(0.25)
percentile_25_balanced_threshold = data.iloc[:, 2].quantile(0.25)
percentile_75_bio_threshold = data.iloc[:, 1].quantile(0.75)
percentile_75_balanced_threshold = data.iloc[:, 2].quantile(0.75)
percentile_97_5_bio_threshold = data.iloc[:, 1].quantile(0.975)
percentile_97_5_balanced_threshold = data.iloc[:, 2].quantile(0.975)

# Save median thresholds overall
out_file_overall=os.path.join(out_path,'FILTERED_median_threshold_per_split.csv')
result_df = pd.DataFrame({'split_type': ['biological', 'balanced'],
                        'median_tpm': [median_bio_threshold, median_balanced_threshold],
                        '2_5percentile_tpm': [percentile_2_5_bio_threshold, percentile_2_5_balanced_threshold],
                        '25percentile_tpm': [percentile_25_bio_threshold, percentile_25_balanced_threshold],
                        '75percentile_tpm': [percentile_75_bio_threshold, percentile_75_balanced_threshold],
                        '97_5percentile_tpm': [percentile_97_5_bio_threshold, percentile_97_5_balanced_threshold],
                        })
result_df.to_csv(out_file_overall, mode='a', header=(not os.path.exists(out_file_overall)), index=False)
