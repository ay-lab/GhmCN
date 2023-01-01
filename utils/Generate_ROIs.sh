#!/usr/bin/env bash
# PreProcess:
# Usage:
#       ./Generate_ROIs.sh  <CellType>
# $1 = cell code name as listed in rnaseq.csv and in the data folder
[[ "$1" == "" ]] && { echo "Two arguments required: <Hi-C window size Bases: default 10000> <CellType: as appears under folder src/data>"; exit 0 ; }

cell_type=$1 # Has to match what is written under the `srt/data/<CellType>`
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUT_DIR=${SCRIPT_DIR%*/}/data/${cell_type}
tmp=`mktemp`
bp=10000

# Check if output folder is present and has content
[[ -d $OUT_DIR ]] || { printf 'ERROR: The specified CellType "<%s>" is not present under the data folder %s\nRun Rice_C.sh or generate appropriate input files\n' ${cell_type} ${SCRIPT_DIR%*/}/data ;  exit 0 ; }
[[ -f $OUT_DIR/hic_chr1.txt ]] ||  { printf 'ERROR: The specified CellType`s data folder "<%s>" has no hic_chr[1-N].txt information\nRun Rice_C.sh or generate appropriate input files\n' $OUT_DIR ;  exit 0 ; }
[[ -f $OUT_DIR/hic_chr18.txt ]] ||  { printf 'ERROR: The specified CellType`s data folder "<%s>" has no hic_chr18.txt file\nIf this is mice, Run Rice_C.sh or generate appropriate input files\nOtherwise modify this file to adjust for your studied genome' $OUT_DIR ;  exit 0 ; }

printf '#!/usr/bin/env R\nlibrary(GenomicRanges)\nlibrary(GenomicAlignments)\noptions(scipen=999)\nROI = rbind(\n' > ${tmp}_Generate_ROIs_mm10_${cell}_${bp}bp.R
# Since we used mouse genome, the chromosomes we consider are up to the 18th
# mod this for loop and the lines after to adjust to your organism (21th for human)
for i in {1..18}; do # <<< MOD this if not mice
 Start=`head $OUT_DIR/hic_chr${i}.txt | cut -f1-2 | tr '\t' '\n' | sort -n | head -n1 `
 End=`tail   $OUT_DIR/hic_chr${i}.txt | cut -f1-2 | tr '\t' '\n' | sort -n | tail -n1 `
 printf 'cbind("chr%s",seq(%s,%s-%s,%s),seq(%s+%s,%s,%s)),\n' $i $Start $End $bp $bp $Start $bp $End $bp >> ${tmp}_Generate_ROIs_mm10_${cell}_${bp}bp.R
done
i=19 # <<< MOD this if not mice
Start=`head $OUT_DIR/hic_chr${i}.txt | cut -f1-2 | tr '\t' '\n' | sort -n | head -n1 `
End=`  tail $OUT_DIR/hic_chr${i}.txt | cut -f1-2 | tr '\t' '\n' | sort -n | tail -n1 `
printf 'cbind("chr%s",seq(%s,%s-%s,%s),seq(%s+%s,%s,%s)))\n
counts = ROI
ROI = GRanges(seqnames=ROI[,1], ranges=IRanges(start=as.numeric(ROI[,2]), end=as.numeric(ROI[,3])))\n
save(list=c("counts","ROI"),file="%s/ROI_mm10_%s_%sbp.RData")\n' $i $Start $End $bp $bp $Start $bp $End $bp $OUT_DIR $cell $bp  >> ${tmp}_Generate_ROIs_mm10_${cell}_${bp}bp.R

R CMD BATCH ${tmp}_Generate_ROIs_mm10_${cell}_${bp}bp.R ${tmp}_Generate_ROIs_mm10_${cell}_${bp}bp.Rout
# cat $Rfolder/Generate_ROIs_mm10_${cell}_${bp}bp.Rout
