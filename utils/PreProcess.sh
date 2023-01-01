#!/usr/bin/env bash
# PreProcess:
# Usage:
#       ./PreProcess.sh  <DIR> <CELL> <TotalEpiMarks> <ChIPBINsize> 
# $1 = folder path
# $2 = cell code name as listed in rnaseq.csv 
# $3 = How many EpiMarks are bein gintegrated? 1, 2, ...,
# $4 = Window size (paper used 10000)
# $5 = Regression flag (1=regression, 0=classification task)
[[ "$1" == "" ]] && { echo "Five arguments required: <DIR> <CELL> <TotalEpiMarks> <ChIPBINsize> <RegressionFlag>"; exit 0 ; }
[[ "$2" == "" ]] && { echo "Five arguments required: <DIR> <CELL> <TotalEpiMarks> <ChIPBINsize> <RegressionFlag>"; exit 0 ; }
[[ "$3" == "" ]] && { echo "Five arguments required: <DIR> <CELL> <TotalEpiMarks> <ChIPBINsize> <RegressionFlag>"; exit 0 ; }
[[ "$4" == "" ]] && { echo "Five arguments required: <DIR> <CELL> <TotalEpiMarks> <ChIPBINsize> <RegressionFlag>"; exit 0 ; }
[[ "$5" == "" ]] && { echo "Five arguments required: <DIR> <CELL> <TotalEpiMarks> <ChIPBINsize> <RegressionFlag>"; exit 0 ; }
[[ -f $1/PreProcess_rf1.out ]] && rm -rf $1/PreProcess_rf1.out 
[[ -f $1/PreProcess_rf0.out ]] && rm -rf $1/PreProcess_rf0.out 

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

conda activate GhmCN
cd ${SCRIPT_DIR%/*}/src
python process_inputs_EGA.py -c $2 -hm $3 -cr $4 -rf $5  -hc 21 
python run_models_EGA.py     -c $2 -hm $3 -cr $4 -rf $5         -df data 
