#!/usr/bin/env R

args = commandArgs(trailingOnly=TRUE)
RNA = read.table(args[1],header=TRUE) 
RNA = RNA[,c(1,6:7)]
RPK = RNA[,3]/(RNA[,2]/1000)
RPM = sum(RPK)/1000000
TPM = round(RPK/RPM,3)
write.table(TPM,file=args[2],row.names=FALSE,col.names=FALSE,quote=FALSE)
