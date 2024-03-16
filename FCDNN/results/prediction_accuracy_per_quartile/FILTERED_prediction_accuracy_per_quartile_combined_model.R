# NOTE:
# File with the sample distributions
# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN//data/sample_groups.csv
# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/prediction_accuracy_per_quartile; R
library(ggplot2)
library(viridis)
library(dplyr)

options(scipen = 999, width = 240)
output_path = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/prediction_accuracy_per_quartile'
data_path = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/'
all_samples_tpm = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences/FILTERED_all_samples_tpm.csv", sep=",", header = TRUE)
quartiles_per_sample <- apply(all_samples_tpm[,-1], 2, quantile, probs = c( 0.25, 0.5, 0.75))

accuracy_per_quartile <- function(df, model, total_genes=1340, quartiles_per_sample){
  quartiles=c("Q1","Q2","Q3","Q4")
  Results = matrix(,ncol=3,nrow=0)
  # Add quartile
  df$Quartile = "Q4"
  # Iterate through blocks of data (samples had 1340 genes)
  for(i in 1:(dim(df)[1]/total_genes)){
    df$Quartile = "Q4"
    range_i = (total_genes*(i-1)+1):(total_genes*i)
    current_quantiles = as.numeric(c(quartiles_per_sample[,i]))
    # Mark Quartiles
    df$Quartile[range_i][df$TPM[range_i] <= current_quantiles[3] ] = "Q3"
    df$Quartile[range_i][df$TPM[range_i] <= current_quantiles[2] ] = "Q2"
    df$Quartile[range_i][df$TPM[range_i] <= current_quantiles[1] ] = "Q1"
    df_temp = df[range_i,]
    # Calc Acc per quartil
    acc=c()
    for(q in quartiles){
      q1 = df_temp$Quartile == q
      acc = c(acc, mean(df_temp[q1,]$Observed == df_temp[q1,]$Predicted))
    }
    # Append to results
    Results = rbind(Results,cbind(acc, quartiles, model))
  }
  return(Results)
}

# Process Combined Model
predictions = read.csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/rerun/C2B/testing/merged/predictions_CleanSamples.csv", header = TRUE)[,1:2]
df = read.table("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_Test_M.txt", header = FALSE)[,c(3,6)]
all(df[,1] == predictions$Observed) # Check that all match
df = cbind(df, predictions$Predicted)[,c(1,3,2)] # Merge and clean
colnames(df)=c("Observed", "Predicted", "TPM") #rename
model="Combined"
results = accuracy_per_quartile(df, model, quartiles_per_sample = quartiles_per_sample)
mean(predictions$Observed == predictions$Predicted)
# Make Plot
Results = data.frame(
  Accuracy = as.numeric(results[,1]),
  Quartile = as.factor(results[,2]),
  Model = as.factor(results[,3]) 
  )

  group.colors = c("Combined" = "black")
  pdf(paste0(output_path,"/FILTERED_prediction_accuracy_per_quartile_combined_model.pdf"),width=4)
    ggplot(Results, aes(Quartile, Accuracy, color=Model)) +
      geom_boxplot(size=1.4) +
      geom_jitter(width = 0.1, alpha = 0.5) +  # Add jittered points with some transparency
      scale_color_manual(values=group.colors) +
      ggtitle("Prediction Accuracy per Quartile") +
      ylab("Accuracy distribution") +
      ylim(0,1) +
      xlab("") +
    theme_bw() +
    theme_minimal()
  dev.off()
