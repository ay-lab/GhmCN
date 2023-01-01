#!/usr/bin/env bash

# Script: Rice_C.sh
# Reformat ice'd normalized Hi-C matrices.
# ./Rice_C.sh <ICED matrix file> <RAW Reference bed> <CellType>
[[ -f $1 ]] || { printf '     [ICED] file %s not found\n' $1 ; exit 0 ; }
[[ -f $2 ]] || { printf '[REFERENCE] file %s not found\n' $2 ; exit 0 ; }

# Make the outdir self-contained in the appropriate fodler under repo/src/data
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUT_DIR=${SCRIPT_DIR%*/}/data/${3}
[[ -d $OUT_DIR ]] && || { printf 'WARNING: the specified cellType`s folder already exists\n' ; } || mkdir $OUT_DIR

tmp=`mktemp`
printf 'TMP file: %s\nSplitting into columns\n' $tmp
split --lines=1 --suffix-length=1  $1 $tmp
tr ' ' '\n' < ${tmp}a > ${tmp}_Line1.txt &
tr ' ' '\n' < ${tmp}b > ${tmp}_Line2.txt &
tr ' ' '\n' < ${tmp}c > ${tmp}_Line3.txt &
wait

printf 'Splitting into lines\n'
split --lines=1000000 --suffix-length=4 --numeric-suffixes ${tmp}_Line1.txt ${tmp}_SpLine1_ &
split --lines=1000000 --suffix-length=4 --numeric-suffixes ${tmp}_Line2.txt ${tmp}_SpLine2_ &
split --lines=1000000 --suffix-length=4 --numeric-suffixes ${tmp}_Line3.txt ${tmp}_SpLine3_ &
wait

printf 'R transform\n'
totFiles=`ls ${tmp}_SpLine3_* | tail -n1 | rev | cut -f1 -d _ | rev `
printf '
#!/usr/bin/env Rscript
library(stringr)
options(scipen=999, width=230)
for(x in 0:%s){
 if(!x%%%%50){cat(x,"/",%s,"\\n")}
 for(y in 1:3){
  A = read.table(       paste0("%s",y,"_",str_pad(x, 4, pad = "0")),header=FALSE,stringsAsFactors=FALSE)
  write.table(A, file = paste0("%s",y,"_",str_pad(x, 4, pad = "0"),".txt"), row.names=FALSE,col.names=FALSE,quote=FALSE)
 }
}' $totFiles $totFiles ${tmp}_SpLine ${tmp}_SpLine > ${tmp}.R
R CMD BATCH ${tmp}.R ${tmp}.out

printf 'Catenating columns\n'
cat ${tmp}_SpLine1_[0-9][0-9][0-9][0-9].txt > ${tmp}_Column01.txt &
cat ${tmp}_SpLine2_[0-9][0-9][0-9][0-9].txt > ${tmp}_Column02.txt &
cat ${tmp}_SpLine3_[0-9][0-9][0-9][0-9].txt > ${tmp}_Column03.txt &
wait

printf 'Pasting columns\n'
paste ${tmp}_Column0[123].txt > ${tmp}.ReFormat

printf 'HiC_Pro2Readable.pl\n'
perl ${SCRIPT_DIR}/HiC_Pro2Readable.pl -r $2 -q ${tmp}.ReFormat -o ${tmp}_${3}_10000_iced.txt

printf 'Separate by chromosome\n'
awk -v OFS="\t" '{if($1!=$3){next} print $2,$4,$5 > "'$OUT_DIR'/hic_"$1".txt"}' ${tmp}_${3}_10000_iced.txt
