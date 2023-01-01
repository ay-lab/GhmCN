#!/usr/bin/env bash
# Usage: 
# $Bfolder/AddCellExpression.sh  <CELL_NAME> <FEATURE_COUNTS_FILE>
[[ "$2" == "" ]] && { echo "Two arguments required: <CELL_NAME> <FEATURE_COUNTS_FILE> "; exit 0 ; }
[[ "$3" != "" ]] && { echo "Two arguments required: <CELL_NAME> <FEATURE_COUNTS_FILE> "; exit 0 ; }
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
RNA=${SCRIPT_DIR%/*}/src/data/rnaseq.csv
t=`mktemp`

# Check if present by using ReGex's + sed + grep
IsPresent=`sed '1!d' $RNA | grep $1'.\?' | wc -c`
[[ "$IsPresent" != "0" ]] && { printf "\n\tCell Line '%s' already in registry... Please add a new one of change given name.\n\n" $1 ; exit 0 ; }

# TPM normalize
Rscript ${SCRIPT_DIR}/FeatureCounts2TPM.R $2 $t

# Then, add cell-name on top
sed -i '1 i\'$1'' $t

# Finally, merge with current registry:
paste $RNA $t -d , > ${t}2
mv ${t}2 $RNA

