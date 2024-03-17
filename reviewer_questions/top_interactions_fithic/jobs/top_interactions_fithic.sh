# # NOTES:
# top 10 interaction has 20 thousand edges
# Create another network from fithic interactions that has at most 20 thousand edges.
# - Filter to nodes having genes

# - Active B cells
# - What is the longest interaction that we have in our dataset, or the 99th percentile
# - Make a histogram of those distances


# srun --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=100g --time=06:00:00 --pty bash -i
# Sconda conda activate ghmc1; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/top_interactions_fithic
out_dir=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/top_interactions_fithic
b72_path=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/01.HiC_CasellasPublicDatasets/02.Analysis/MergedData/Bcell_72h/hic_results/matrix/Bcell_72h
fithic_utils=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/fithic/fithic/utils

# Inputs for conversion:
valid_pairs=$b72_path/../../data/Bcell_72h/Bcell_72h.allValidPairs
matrix=$b72_path/raw/10000/Bcell_72h_10000.matrix
bed=$b72_path/raw/10000/Bcell_72h_10000_abs.bed
bias=$b72_path/iced/10000/Bcell_72h_10000_iced.matrix.biases
bias_corrected=$b72_path/iced/10000/Bcell_72h_10000_iced.matrix_corrected.biases
iced_matrix=$b72_path/iced/10000/Bcell_72h_10000_iced.matrix

# # >>>>>> NOTE about the incorrect bias:
# # 
# # When running HiCPro2FitHiC.py for the first time I had an id mismatch error between the bias trying to assing/access info out of range for the matrix data.
# # The node IDs in the matrix data go up to 273121
# grep 273122 $matrix #empty
# head $matrix
# # However, the bias file has 273122 lines
# wc -l $bias
# # Unless the indexing begins at 0, the bias has apparently 1 extra line making the files to not match in numbers of lines and total node IDs
# # 
# # Checking at the bias, the first value occured at line 302
# head $bias
# tail $bias
# head -n302 $fbias | tail
# # and the matrix data starts having information at the node ID 301
# head $matrix
# # Then, the bias line 302 annotates to the node ID 301,
# # But it is being anotated to the node ID 302 instead, reason why I got an error (I believe, have no clear evidence).
# # Thus, I shall remove the first line such that the line 302 becomes the 301
tail -n 273121 $bias > $bias_corrected

# Transform HiCPro to fithic data inputs
mkdir -p $out_dir/hicpro_to_fithic/activ72_B
python $fithic_utils/HiCPro2FitHiC.py \
  -i $matrix \
  -b $bed \
  -s $bias_corrected \
  -r 10000 \
  -o $out_dir/hicpro_to_fithic/activ72_B
# NOTE: When running HiCPro2FitHiC.py for the first time I had an id mismatch between the bias trying to assing/access info out of range for the matrix data. The node IDs in the matrix data go up to 273121. However, the bias file has 273122 lines. Unless the indexing begins at 0, they dont match in numbers of IDs and lines... Checking at the bias, the first value occurred at line 302, and the matrix data starts having information at the node ID 301. Then, the bias line 302 annotates to the node ID 301. But it is being annotated to the node ID 302 instead, the reason why I got an error. Thus, I shall remove the first line such that the line 302 becomes the 301, and the program will (hopefully) assign na bias if I missed one
# e85722ec6ab27b06c27e43b3c65c604d  fithic.biases.gz
# 3595e8b189e4e021937ca657cbcb3153  fithic.fragmentMappability.gz
# ee600477b387d28a6d2e81a27b4cae89  fithic.interactionCounts.gz

# Bash version while we try to make the other work
# Note: This resulted in a corrupted gzip `fithic.interactionCounts.gz` file...
mkdir -p $out_dir/hicpro2fithic/activ72_B
bash $fithic_utils/validPairs2FitHiC-fixedSize.sh 10000 activ72h $valid_pairs $out_dir/hicpro2fithic/activ72_B
# a0101e17f5ef3e30b05b815dd66b3721  activ72h_fithic.contactCounts.gz
# 1339b7c6a098ba1bde925de4a9264f2e  fithic.biases.gz
# 710c2c4d8f1723e75cb714441918ae4c  fithic.fragmentMappability.gz
# fc61e23ee650ad78fc54e2e850fd4494  fithic.interactionCounts.gz



ll $out_dir/hicpro_to_fithic/activ72_B
bias_file=$out_dir/hicpro_to_fithic/activ72_B/fithic.biases.gz
fragment_file=$out_dir/hicpro_to_fithic/activ72_B/fithic.fragmentMappability.gz
interactions_file=$out_dir/hicpro_to_fithic/activ72_B/fithic.interactionCounts.gz
# interactions_file/hicpro_to_fithic/activ72_B/activ72h_fithic.contactCounts.gz
out_fithic=$out_dir/fithic_results/hicpro_to_fithic
mkdir -p $out_fithic
