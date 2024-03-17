# Sconda conda activate TensorFlow_CPU_01; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/deeplift_analysis; R
workdir = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/deeplift_analysis/'
file_data = paste0(workdir,"deeplift_scores.RData")
plot_file = paste0(workdir,"FILTERED_deeplift_clean_samples_scores.pdf")
plot_file_nogrid = paste0(workdir,"FILTERED_deeplift_clean_samples_scores_nogrid.pdf")

library(ggplot2)
library(viridis)
library(plyr)
library(dplyr)
options(width=240,scipen=999)
library(stringr)

# >>> Plot aggregate:
data_summary <- function(data, groupnames){
  var = "Score"
      require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  var)
  data_sum <- rename(data_sum, c(Score = "mean"))
 return(data_sum)
}
if(!file.exists(file_data)){
  # >>> Read deepLift scores
  ZR_TP_B_Net01 = as.matrix(read.csv("TPscores_NeutralRef_CleanSamples_Net200-100-50-1.csv",header=FALSE))
  ZR_TN_B_Net01 = as.matrix(read.csv("TNscores_NeutralRef_CleanSamples_Net200-100-50-1.csv",header=FALSE))

  # >>> Generate names per bin
  Bin = c(
    paste0("Prom_Neg_", str_pad(50:1 , 3, pad = "0")),
    paste0("Prom_Pos_", str_pad(1:50 , 3, pad = "0")),
    paste0("TSS_", str_pad(15:1 , 3, pad = "0")),
    paste0("Gbody_", str_pad(1:100, 3, pad = "0")),
    paste0("TTS_", str_pad(1:15 , 3, pad = "0"))
    )
  
  # >>> Merge scores
  zTP  = matrix(c(ZR_TP_B_Net01), ncol=1, byrow=TRUE)
  bin  = rep(Bin,each=dim(ZR_TP_B_Net01)[1])
  lab  = rep("Neutral_TP", each=dim(ZR_TP_B_Net01)[1]*dim(ZR_TP_B_Net01)[2])
  zTP  = cbind(zTP,bin,lab)

  zTN  = matrix(c(ZR_TN_B_Net01), ncol=1, byrow=TRUE)
  bin  = rep(Bin,each=dim(ZR_TN_B_Net01)[1])
  lab  = rep("Neutral_TN", each=dim(ZR_TN_B_Net01)[1]*dim(ZR_TN_B_Net01)[2])
  zTN  = cbind(zTN,bin,lab)

  All  = rbind(zTP,zTN)

  # Combine into a data frame and reorg naming
  All = data.frame(Score = as.numeric(All[,1]) , Bin = as.factor(All[,2]), Label = as.factor(All[,3]))
  All$Bin   = factor(All$Bin   , levels = levels(All$Bin  )[c(150:101,151:200,215:201,1:100,216:230)])
  All$Label = factor(All$Label , levels = levels(All$Label)[c(1,4,3,2)])

  save(All, file=file_data)
}

load(file_data)
groupnames=c("Bin", "Label")
df2 <- data_summary(All, groupnames=c("Bin", "Label"))
head(df2)

pdf(plot_file,width=28,height=7)
ggplot(df2, aes(x=Bin, y=Score, group=Label, color=Label)) +
           geom_hline(yintercept=0, color="black") +
       geom_errorbar(aes(ymin=Score-sd, ymax=Score+sd), width=.5,
                     position=position_dodge(0.05),alpha=0.35) +
       scale_colour_manual(values=c("darkblue","darkred")) +
       ggtitle("Importance Score Distribution per bin: Clean Samples [49 samples]", subtitle = "Neutral Reference [activation value ~ 0.5]") +
       geom_point(size=2.5)+
       ylim(-1.2,2.85)+
       geom_line(alpha=0.25)+
       theme_bw() +
       theme_minimal() +
       theme(axis.ticks.x = element_blank(),
             axis.text.y  = element_text(face = "bold",size = 16),
             legend.text  = element_text(face = "bold",size = 16),
             legend.title = element_text(face = "bold",size = 16),
             title        = element_text(face = "bold",size = 20)) +
       theme(axis.text.x = element_text(angle = 90)) 
dev.off()

pdf(plot_file_nogrid,width=28,height=7)
ggplot(df2, aes(x=Bin, y=Score, group=Label, color=Label)) +
           geom_hline(yintercept=0, color="black") +
       geom_errorbar(aes(ymin=Score-sd, ymax=Score+sd), width=.5,
                     position=position_dodge(0.05),alpha=0.35) +
       scale_colour_manual(values=c("darkblue","darkred")) +
       ggtitle("Importance Score Distribution per bin: Clean Samples [49 samples]", 
               subtitle = "Neutral Reference [activation value ~ 0.5]") +
       geom_point(size=2.5)+
       ylim(-1.2,2.85)+
       geom_line(alpha=0.25)+
       theme_minimal() +
       theme(axis.ticks.x = element_blank(),
             axis.text.y  = element_text(face = "bold",size = 16),
             legend.text  = element_text(face = "bold",size = 16),
             legend.title = element_text(face = "bold",size = 16),
             title        = element_text(face = "bold",size = 20),
             axis.text.x = element_text(angle = 90),
             panel.grid.minor.y = element_blank(), 
             panel.grid.major.x = element_blank(), 
             panel.grid.minor.x = element_blank()  # Remove vertical minor grid lines if needed
) 
dev.off()
