# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/traditional_ml_methods; R
#!/usr/bin/R
# All.vs.All.R
wkdir="/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/traditional_ml_methods/"
setwd(wkdir)
library(ggplot2)
library(viridis)
library(reshape2)

if(!file.exists(paste0(wkdir,"traditional_ml_methods_data_formatted.RData"))){
  # Load raw data and reformat to save as melted and ordered
  Freq = read.csv(paste0(wkdir,"FILTERED_auc_per_method_per_sample.csv"))
  Freq <- melt(Freq, id.vars = "Sample", variable.name = "Model", value.name = "AUC")
  # Reorder Model's data
  Freq$Model = factor(Freq$Model, levels = c("DNN","LRg","RFo","SVM"))
  # Reorder Samples by Ascending AUC values
  subO = Freq$Model == "DNN"
  Freq$Order = factor(rep(order(order(Freq$AUC[subO])),4)  , levels=1:table(Freq$Model)[1])
  Freq$Sample = factor(Freq$Sample, levels = Freq$Sample[subO][order(Freq$AUC[subO])])
  # Declare specific color used
  group.colors = c(SVM = "darkgrey" ,RFo = "dimgrey" ,LRg = "gold4" ,DNN = "darkgreen")
  # Save data as table and RData
  write.table(file = paste0(wkdir,"traditional_ml_methods_data_formatted.csv"), x=Freq, quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE)

  save(Freq, group.colors, file = paste0(wkdir,"traditional_ml_methods_data_formatted.RData"))
}
load(file = paste0(wkdir,"traditional_ml_methods_data_formatted.RData"))

pdf(paste0(wkdir,"FILTERED_traditional_ml_balanced.pdf"),height=5,width=4)
  ggplot(Freq, aes( Model, AUC, color=Model)) +
    geom_boxplot(size=1.4) +
    geom_jitter(width = 0.1, alpha = 0.5) +  # Add jittered points with some transparency
    scale_color_manual(values=group.colors) +
    ggtitle("") +
    ylab(" ") +
    xlab("") +
  theme_bw() +
  theme_minimal()
dev.off()

pdf(paste0(wkdir,"FILTERED_traditional_ml_balanced_bar.pdf"),height=10, width=4)
  ggplot(Freq, aes( Order, AUC, fill=Model)) +
    geom_bar(linewidth=0.05,stat="identity", position=position_dodge()) +
    scale_fill_manual(values=group.colors) +
    ggtitle("AUC Score Distribution per Model\n(Balanced)") +
    ylab("AUC Score") +
    xlab("") +
    ylim(0,1) +
  theme_bw() +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

pdf(paste0(wkdir,"FILTERED_traditional_ml_balanced_bar.pdf"),height=10, width=8)
  ggplot(Freq, aes( Sample, AUC, fill=Model)) +
    geom_bar(linewidth=0.05,stat="identity", position=position_dodge()) +
    scale_fill_manual(values=group.colors) +
    ggtitle("AUC Score Distribution per Sample\n(Balanced)") +
    ylab("AUC Score") +
    xlab("") +
    ylim(0,1) +
  theme_bw() +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=5),
        axis.ticks.y=element_blank())
dev.off()

q('no')

# Test for normal distribution
for(model in levels(Freq$Model)){
  if (shapiro.test(Freq[Freq$Model==model,]$AUC)$p.value < 0.05){print(paste0(model," not normal"))}
}
# Print the result
model = "DNN"
shapiro.test(Freq[Freq$Model==model,]$AUC)
model = "RFo"
shapiro.test(Freq[Freq$Model==model,]$AUC)

# Test for normal distribution
for(model in levels(Freq$Model)[-1]){
    test_result <- wilcox.test(Freq[Freq$Model=="DNN",]$AUC, Freq[Freq$Model==model,]$AUC, paired = TRUE)
    print(paste0("DNN vs ",model," p-value: ",test_result$p.value))
}
# Print the result
model = "DNN"
shapiro.test(Freq[Freq$Model==model,]$AUC)
model = "RFo"
shapiro.test(Freq[Freq$Model==model,]$AUC)
model = "LRg"
shapiro.test(Freq[Freq$Model==model,]$AUC)
model = "SVM"
shapiro.test(Freq[Freq$Model==model,]$AUC)