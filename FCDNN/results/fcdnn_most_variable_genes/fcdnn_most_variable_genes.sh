# This file contains on the 4th column the most variable genes.
#   Use it to go into the metadata of the CleanSamples model and get where these genes are locatet at from Train/Test/Dev datasets
variable_genes_file=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences/FILTERED_variable_label_genes_promoter_coords.bed

# CleanSamples location:
clean_f=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_
clean_dev=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_Dev_M.txt
clean_test=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_Test_M.txt
clean_train=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_Train_M.txt

cut -f4 $variable_genes_file | sort > /tmp/variable_genes.txt
sort -k1,1 $clean_dev > /tmp/dev.txt

for i in Dev Test Train; do
    f=${clean_f}${i}_M.txt
    join /tmp/variable_genes.txt <(sort -k1,1 $f | cut -f1) | uniq > /tmp/${i}.txt 
    wc -l  /tmp/${i}.txt
done

cat /tmp/[DT]*.txt  > /tmp/all.txt

# 150 /tmp/Dev.txt
# 145 /tmp/Test.txt
# 2152 /tmp/Train.txt

#  All genes found
# I can do a inner merge in python and get the index numbers of the matchs of the M files and then use that to get the needed coluns of each X/Y file
