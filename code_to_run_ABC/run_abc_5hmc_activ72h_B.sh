# To Start a GPU-interactive job:
# srun --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=100g --time=04:00:00 --pty bash -i

name=activ72h_B
rep=rep1
peak_name=${name}_${rep}
g2a_input=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input
g2a_output_orig=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_output
g2a_output=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_output_5hmc
atac_bam=$g2a_input/${name}_atac_rep1.bam
atac_bam_rep2=$g2a_input/${name}_atac_rep2.bam
k27ac_bam=$g2a_input/${name}_h3k27ac.bam
hmc_bam=$g2a_input/${name}_5hmc_r1.bam
hmc_bam_r2=$g2a_input/${name}_5hmc_r3.bam

### Step 1. Define candidate elemets
# Taken from prev: same ATAC files

### Step 2. Quantifying Enhancer Activity:
# Sconda
cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ABC
conda deactivate; conda activate abc

python src/run.neighborhoods.py \
--candidate_enhancer_regions $g2a_output_orig/Peaks/${peak_name}.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes $g2a_input/genes_tss_tts.bed \
--H3K27ac ${hmc_bam},${hmc_bam_r2} \
--DHS ${atac_bam},${atac_bam_rep2} \
--expression_table $g2a_input/${name}_expression.tsv \
--chrom_sizes $g2a_input/mm10.chrom.sizes.extra \
--ubiquitously_expressed_genes $g2a_input/ubiquitous_expressed_genes.csv \
--cellType ${name} \
--outdir $g2a_output/Neighborhoods/ 

# With h3k27ac
# Total enhancers: 149547
#          Promoters: 23623
#          Genic: 64039
#          Intergenic: 61885
# 
# With 5hmC
# Total enhancers: 149547
#          Promoters: 23623
#          Genic: 64039
#          Intergenic: 61885

### Step 3. Computing the ABC Score
# Compute ABC scores by combining Activity (as calculated by `run.neighborhoods.py`) and Hi-C.
cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ABC
python src/predict.py \
  --enhancers $g2a_output/Neighborhoods/EnhancerList.txt \
  --genes $g2a_output/Neighborhoods/GeneList.txt \
  --HiCdir $g2a_input/hic/${name}/iced \
  --hic_type bedpe \
  --chrom_sizes $g2a_input/mm10.chrom.sizes \
  --hic_resolution 10000 \
  --scale_hic_using_powerlaw \
  --threshold .02 \
  --cellType ${name} \
  --outdir $g2a_output/Predictions/ \
  --make_all_putative
