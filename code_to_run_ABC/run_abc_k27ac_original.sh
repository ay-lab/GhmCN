# To Start a GPU-interactive job:
# srun --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=100g --time=04:00:00 --pty bash -i

name=activ72h_B
rep=rep1
peak_name=${name}_${rep}
g2a_input=/mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/ghmc_to_abc/ABC_input
g2a_output=/mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/ghmc_to_abc/ABC_output
atac_bam=$g2a_input/${name}_atac_rep1.bam
atac_bam_rep2=$g2a_input/${name}_atac_rep2.bam
k27ac_bam=$g2a_input/${name}_h3k27ac.bam

### Step 1. Define candidate elemets
Sconda
conda deactivate; conda activate macs
#Call peaks
macs2 callpeak \
-t ${atac_bam} \
-n ${peak_name}.macs2 \
-f BAM \
-g mm \
-p .1 \
--call-summits \
--outdir $g2a_output/Peaks/ 

#Sort narrowPeak file
bedtools sort -faidx $g2a_input/mm10.chrom.sizes.extra -i $g2a_output/Peaks/${peak_name}.macs2_peaks.narrowPeak > $g2a_output/Peaks/${peak_name}.macs2_peaks.narrowPeak.sorted
wc -l $g2a_output/Peaks/${peak_name}.macs2_peaks.narrowPeak.sorted
# 3184319 rep1
# 3083088 rep2
# 3366838 new_rep1

#Call candidate regions
# NOTE: missed `pyranges` & `pysam` package
# conda install pyranges pysam
cd /mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/ABC
python src/makeCandidateRegions.py \
--narrowPeak $g2a_output/Peaks/${peak_name}.macs2_peaks.narrowPeak.sorted \
--bam ${atac_bam_rep2} \
--outDir $g2a_output/Peaks/ \
--chrom_sizes $g2a_input/mm10.chrom.sizes.extra \
--regions_blocklist $g2a_input/mm10-blacklist.v2.bed \
--regions_includelist $g2a_input/whitelist_genes_tss.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000 
wc -l $g2a_output/Peaks/${peak_name}.macs2_peaks.narrowPeak.sorted.candidateRegions.bed
# 149547

### Step 2. Quantifying Enhancer Activity:
conda deactivate; conda activate abc

python src/run.neighborhoods.py \
--candidate_enhancer_regions $g2a_output/Peaks/${peak_name}.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes $g2a_input/genes_tss_tts.bed \
--H3K27ac ${k27ac_bam} \
--DHS ${atac_bam},${atac_bam_rep2} \
--expression_table $g2a_input/${name}_expression.tsv \
--chrom_sizes $g2a_input/mm10.chrom.sizes.extra \
--ubiquitously_expressed_genes $g2a_input/ubiquitous_expressed_genes.csv \
--cellType ${name} \
--outdir $g2a_output/Neighborhoods/ 
# rep1
# Total enhancers: 149547
#          Promoters: 23623
#          Genic: 64039
#          Intergenic: 61885

### Step 3. Computing the ABC Score
# Compute ABC scores by combining Activity (as calculated by ```run.neighborhoods.py```) and Hi-C.
cd /mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/ABC
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

### Step 4. Get Prediction Files for Variant Overlap 
- this step is already being performed in src/predict.py

Sample Command:
```
python src/getVariantOverlap.py \
--all_putative example_chr22/ABC_output/Predictions/EnhancerPredictionsAllPutative.txt.gz \
--chrom_sizes example_chr22/reference/chr22 \
--outdir .
```
