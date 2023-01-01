#!/usr/bin/env R

args = commandArgs(trailingOnly=TRUE)
if( length(args)!=5) {   stop("\tThe 4 arguments are required: <BAM> <CONDITION_NAME> <OUTDIR> <CELLNAME>\n\n", call.=FALSE) }

library(GenomicRanges)
library(GenomicAlignments)
options(scipen=999)


bam  = args[1]
cond = args[2]
cell = args[3]
out  = args[4]
bp   = 10000
file2Load = paste0(out,"/ROI_mm10_",cell,"_",bp,"bp.RData")
cat("\n# >>> [",cell,"] Loading ROI and cound objects\n# file: ",file2Load)
load(file2Load)
cat("\n# >>> [",cell,"] Reading ROI Signal from BAM file\n# file: ",bam)
Signal = summarizeOverlaps(ROI, bam,inter.feature = FALSE)
colnames(Signal@assays$data@listData$counts) = cond
counts = cbind( counts , Signal@assays$data@listData$counts[,1] )
cat("\n# >>> [",cell,"] Writting ROI Signal to OUT folder:\n")
for(chr in unique(counts[,1])){
  cat(paste0("  # > [",cell,"] Processing ",chr,"\n"))
  outname=paste0(out,"/",chr,"_",bp,"bp_",cond,".count")
  write.table(counts[(counts[,1]==chr),],file = outname,row.names=FALSE,col.names=FALSE,quote=FALSE, sep="\t")
}
