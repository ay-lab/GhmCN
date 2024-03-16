#!/bin/bash
#SBATCH -J fithic_ice_bias_test
#SBATCH -o /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/top_interactions_fithic/run_fithic_ice_bias_test_bias.out
#SBATCH -e /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/top_interactions_fithic/run_fithic_ice_bias_test_bias.err
#SBATCH -t 12:00:00
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 8

source /mnt/BioAdHoc/Groups/RaoLab/Edahi/00.Scripts/01.Downloaded/anaconda3/etc/profile.d/conda.sh
conda activate ghmc1

out_dir=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/top_interactions_fithic

interactions_file=$out_dir/hicpro_to_fithic/activ72_B/fithic.interactionCounts.gz
bias_file=$out_dir/hicpro_to_fithic/activ72_B/fithic.biases.gz
fragment_file=$out_dir/hicpro_to_fithic/activ72_B/fithic.fragmentMappability.gz
out_fithic=$out_dir/fithic_results/hicpro_to_fithic_test_bias
mkdir -p $out_fithic
cd $out_fithic

fithic \
-i $interactions_file \
-t $bias_file \
-f $fragment_file \
-o $out_fithic \
-U 5000000 \
-L 10000 \
-v \
--biasLowerBound 0.25 \
--biasUpperBound 3 \
-r 10000
