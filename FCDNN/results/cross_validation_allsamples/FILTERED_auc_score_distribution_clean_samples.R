# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples; R
#!/usr/bin/R
options(width=240)
library(ggplot2)
library(viridis)

wkdir = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/'
setwd(wkdir)

group.colors  = c(Each_Sample = "black", Cross_Val = "red")

# 1st plot
# NOTES:
# - This was considered an infair comparision because we have 19 chromosomes in the cross-validation
# - whereas we have 49 in the chean Samples model
Freq = read.csv(paste0(wkdir,"FILTERED_auc_score_distribution_clean_samples_cross_val.csv"))
pdf(paste0(wkdir,"FILTERED_clean_samples_model_in_cell_specific_data_and_cross_validation_auc.pdf"),width=6)
  ggplot(Freq, aes( Tested_Dataset, AUCscore, color=Tested_Dataset)) +
    geom_boxplot(size=1.1) +
    geom_jitter(width = 0.1, alpha = 0.5) +  # Add jittered points with some transparency
    scale_color_manual(values=group.colors) +
    ggtitle("AUC Score Distrib. of 'CleanSamples' Model") +
    ylab("AUC Score") +
    xlab("") +
  theme_bw() +
  theme_minimal()
dev.off()

# Only the clean data model
Freq = as.data.frame(Freq[Freq$Tested_Dataset != "Cross_Val",])
pdf(paste0(wkdir,"FILTERED_clean_samples_model_in_cell_specific_data_auc.pdf"),width=3)
  ggplot(Freq, aes( Tested_Dataset, AUCscore, color=Tested_Dataset)) +
    geom_boxplot(size=1.1) +
    geom_jitter(width = 0.1, alpha = 0.5) +  # Add jittered points with some transparency
    scale_color_manual(values=group.colors) +
    ggtitle("AU ROC Score Distrib\nCombined Model") +
    ylab("AUC Score") +
    xlab("") +
  theme_bw() +
  theme_minimal()
dev.off()
# 2nd plot
# - Will plot the AUC distribution per chromosome showing the 49 dots per boxplot, add a line with CleanSamples's median AUC
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/FILTERED_auc_score_distribution_clean_samples.py

# 3rd plot
# - Will plot the AUC distribution per sample per chromosome, that way I can overlay the Clean Samples' sample AUC 
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cross_validation_allsamples/FILTERED_auc_score_distribution_clean_samples.py