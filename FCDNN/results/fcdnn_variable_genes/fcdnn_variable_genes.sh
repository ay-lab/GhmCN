# This file contains on the 4th column the most variable genes.
#   Use it to go into the metadata of the CleanSamples model and get where these genes are locatet at from Train/Test/Dev datasets
variable_genes_file=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences/FILTERED_any_degree_of_variable_label_genes_promoter_coords.bed

# CleanSamples location:
clean_f=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_
clean_dev=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_Dev_M.txt
clean_test=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_Test_M.txt
clean_train=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_Train_M.txt

cut -f4 $variable_genes_file | sort > /tmp/any_variable_genes.txt

for i in Dev Test Train; do
    f=${clean_f}${i}_M.txt
    join /tmp/any_variable_genes.txt <(sort -k1,1 $f | cut -f1) | uniq > /tmp/any_variable_genes_${i}.txt 
    wc -l  /tmp/any_variable_genes_${i}.txt 
done

cat /tmp/any_variable_genes_[DT]*.txt  > /tmp/all_any_variable_genes.txt
wc -l /tmp/all_any_variable_genes.txt
# 658 /tmp/any_variable_genes_Dev.txt
# 652 /tmp/any_variable_genes_Test.txt
# 9752 /tmp/any_variable_genes_Train.txt

#  All genes found
# I can do a inner merge in python and get the index numbers of the matchs of the M files and then use that to get the needed coluns of each X/Y file
