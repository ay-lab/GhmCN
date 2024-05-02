name=activ72h_B
name_alt=resting_B
g2a=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc
g2a_input=$g2a/ABC_input
src=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/src
g2a_hic=$g2a/ABC_input/hic/${name}
g2a_hic_alt=$g2a/ABC_input/hic/${name_alt}
mkdir -p $g2a_hic $g2a_hic_alt/raw $g2a_hic_alt/iced
# Objective:
# 1
# DONE
# Get a BED-ish file that defines the TSS-TTS for a gene and provides its orientation.
# chrom start end name orientation
awk -v OFS="\t" '{t=$4; $4=$5"\t"0; $5=t; print}' /mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/mm10/Genes.gtf | grep -v chrUn_JH584304 > $g2a_input/genes_tss_tts.bed
wc -l $g2a_input/genes_tss_tts.bed

# 2
# DONE
# Also, download Kundaje's blaclisted regions for mm10 to exclude such regions.
# https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
wc -l $g2a_input/mm10-blacklist.v2.bed

# 3
# DONE
# Record paths to B72's 5hmC, ATACseq and H3K27ac's BAM files
cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/03.5hmC_as_Enhancer_Marker/02.Mapping/aB_wt_H3K27Ac/[ar]B_wt_H3K27Ac.mapped.sorted.bam*  $g2a_input/
cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/Bcells/CMSIP_PMID_31028100_Jerry_2019_ScienceImmunology/02.Mapping/WT.72.IP.E1/WT.72.IP.E1.mm10.bam $g2a_input/${name}_5hmc_rep1.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/Bcells/CMSIP_PMID_31028100_Jerry_2019_ScienceImmunology/02.Mapping/WT.72.IP.E4/WT.72.IP.E4.mm10.bam $g2a_input/${name}_5hmc_rep2.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Vipul_and_Jerry/01.CSR_project/05.CMSIP/04.bsmap_deleating_mapped_reads_puro_lambda/WT-72-IP-E1_sorted.bam.bam $g2a_input/${name}_5hmc_r1.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Vipul_and_Jerry/01.CSR_project/05.CMSIP/04.bsmap_deleating_mapped_reads_puro_lambda/WT-72-IP-E2_sorted.bam.bam $g2a_input/${name}_5hmc_r2.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Vipul_and_Jerry/01.CSR_project/05.CMSIP/04.bsmap_deleating_mapped_reads_puro_lambda/WT-72-IP-E3_sorted.bam.bam $g2a_input/${name}_5hmc_r3.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Vipul_and_Jerry/01.CSR_project/05.CMSIP/04.bsmap_deleating_mapped_reads_puro_lambda/WT-0-IP-E1_sorted.bam.bam $g2a_input/${name_alt}_5hmc_r1.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Vipul_and_Jerry/01.CSR_project/05.CMSIP/04.bsmap_deleating_mapped_reads_puro_lambda/WT-0-IP-E2_sorted.bam.bam $g2a_input/${name_alt}_5hmc_r2.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Vipul_and_Jerry/01.CSR_project/05.CMSIP/04.bsmap_deleating_mapped_reads_puro_lambda/WT-0-IP-E3_sorted.bam.bam $g2a_input/${name_alt}_5hmc_r3.bam
# cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/Bcells/ATACSeq_PMID_31028100_Jerry_2019_ScienceImmunology/02.Mapping/WT_[07][02]h_rep2/WT_[07][02]h_rep2.whitelist.bam*  $g2a_input/
# mv $g2a_input/aB_wt_H3K27Ac.mapped.sorted.bam $g2a/${name}_h3k27ac.bam
# mv $g2a_input/rB_wt_H3K27Ac.mapped.sorted.bam $g2a/${name_alt}_h3k27ac.bam
# mv $g2a_input/WT_00h_rep1.whitelist.bam $g2a/${name_alt}_atac_rep1.bam
# mv $g2a_input/WT_00h_rep1.whitelist.bam.bai $g2a/${name_alt}_atac_rep1.bam.bai
# mv $g2a_input/WT_00h_rep2.whitelist.bam $g2a/${name_alt}_atac_rep2.bam
# mv $g2a_input/WT_00h_rep2.whitelist.bam.bai $g2a/${name_alt}_atac_rep2.bam.bai
# cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input/resting_B_atac.bam  $g2a_input/${name_alt}_atac_rep1.bam
# cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input/resting_B_atac.bam.bai  $g2a_input/${name_alt}_atac_rep1.bam.bai
# mv $g2a_input/WT_72h_rep1.whitelist.bam $g2a/${name}_atac_rep1.bam
# mv $g2a_input/WT_72h_rep1.whitelist.bam.bai $g2a/${name}_atac_rep1.bam.bai
# mv $g2a_input/WT_72h_rep2.whitelist.bam $g2a/${name}_atac_rep2.bam
# mv $g2a_input/WT_72h_rep2.whitelist.bam.bai $g2a/${name}_atac_rep2.bam.bai
# Copy fresh H3K4me1 datasets:
cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/mouse_NaiveTcells/ChIPSeq_PMID_31873292_Singh_2019_NatImmunology//02.Mapping/H3K4me1_aB_1/H3K4me1_aB_1.bam $g2a_input/${name}_h3k4me1_r1.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/mouse_NaiveTcells/ChIPSeq_PMID_31873292_Singh_2019_NatImmunology//02.Mapping/H3K4me1_aB_2/H3K4me1_aB_2.bam $g2a_input/${name}_h3k4me1_r2.bam
cp /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/mouse_NaiveTcells/ChIPSeq_PMID_31873292_Singh_2019_NatImmunology//02.Mapping/H3K4me1_rB_1/H3K4me1_rB_1.bam $g2a_input/${name_alt}_h3k4me1_r1.bam


samtools index -@ 8 $g2a_input/${name}_h3k27ac.bam
samtools index -@ 8 $g2a_input/resting_B_h3k27ac.bam
samtools index -@ 8 $g2a_input/${name}_5hmc_rep1.bam
samtools index -@ 8 $g2a_input/${name}_5hmc_rep2.bam
samtools index -@ 8 $g2a_input/${name}_5hmc_r1.bam
samtools index -@ 8 $g2a_input/${name}_5hmc_r2.bam
samtools index -@ 8 $g2a_input/${name}_5hmc_r3.bam
samtools index -@ 24 $g2a_input/${name_alt}_5hmc_r1.bam
samtools index -@ 24 $g2a_input/${name_alt}_5hmc_r2.bam
samtools index -@ 24 $g2a_input/${name_alt}_5hmc_r3.bam
samtools index -@ 24 $g2a_input/${name}_h3k4me1_r1.bam
samtools index -@ 24 $g2a_input/${name}_h3k4me1_r2.bam
samtools index -@ 24 $g2a_input/${name_alt}_h3k4me1_r1.bam


science_data=/mnt/BioAdHoc/Groups/RaoLab/Vipul_and_Jerry/01.CSR_project/01.ATACseq/01.First_Replicate/
# cp ${science_data}/01/02.Mapping/01_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam  $g2a/${name_alt}_atac_rep1.bam
# cp ${science_data}/01/02.Mapping/01_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam.bai $g2a/${name_alt}_atac_rep1.bam.bai
# cp ${science_data}/02/02.Mapping/02_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam  $g2a/${name_alt}_atac_rep2.bam
# cp ${science_data}/02/02.Mapping/02_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam.bai $g2a/${name_alt}_atac_rep2.bam.bai
cp ${science_data}/13/02.Mapping/13_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam   $g2a/${name}_atac_rep1.bam
cp ${science_data}/13/02.Mapping/13_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam.bai  $g2a/${name}_atac_rep1.bam.bai
cp ${science_data}/14/02.Mapping/14_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam   $g2a/${name}_atac_rep2.bam
cp ${science_data}/14/02.Mapping/14_mm10_merge_onlymapped_sorted_rmdup_nochrM.bam.bai  $g2a/${name}_atac_rep2.bam.bai

# Make bw tracks
conda deactivate
conda activate deeptools
# 5hmC bigwig
bamCoverage -b $g2a_input/${name}_5hmc_rep1.bam -o $g2a_input/${name}_5hmc_rep1.bw
bamCoverage -b $g2a_input/${name}_5hmc_rep2.bam -o $g2a_input/${name}_5hmc_rep2.bw
bamCoverage -b $g2a_input/${name}_5hmc_r1.bam -o $g2a_input/${name}_5hmc_r1.bw
bamCoverage -b $g2a_input/${name}_5hmc_r2.bam -o $g2a_input/${name}_5hmc_r2.bw
bamCoverage -b $g2a_input/${name}_5hmc_r3.bam -o $g2a_input/${name}_5hmc_r3.bw
bamCoverage -b $g2a_input/${name_alt}_5hmc_r1.bam -o $g2a_input/${name_alt}_5hmc_r1.bw --numberOfProcessors 24
bamCoverage -b $g2a_input/${name_alt}_5hmc_r2.bam -o $g2a_input/${name_alt}_5hmc_r2.bw --numberOfProcessors 24
bamCoverage -b $g2a_input/${name_alt}_5hmc_r3.bam -o $g2a_input/${name_alt}_5hmc_r3.bw --numberOfProcessors 24
# atac bigwig
bamCoverage -b $g2a_input/${name}_atac_rep1.bam -o $g2a_input/${name}_atac_rep1.bw
bamCoverage -b $g2a_input/${name}_atac_rep2.bam -o $g2a_input/${name}_atac_rep2.bw
bamCoverage -b $g2a_input/${name_alt}_atac_rep1.bam -o $g2a_input/${name_alt}_atac_rep1.bw --numberOfProcessors 24
bamCoverage -b $g2a_input/${name_alt}_atac_rep2.bam -o $g2a_input/${name_alt}_atac_rep2.bw --numberOfProcessors 24
# k27 bigwig
bamCoverage -b $g2a_input/${name}_h3k27ac.bam -o $g2a_input/${name}_h3k27ac.bw
bamCoverage -b $g2a_input/${name_alt}_h3k27ac.bam -o $g2a_input/${name_alt}_h3k27ac.bw  --numberOfProcessors 24
# h3k4me1
bamCoverage -b $g2a_input/${name}_h3k4me1_r1.bam -o $g2a_input/${name}_h3k4me1_r1.bw  --numberOfProcessors 24
bamCoverage -b $g2a_input/${name}_h3k4me1_r2.bam -o $g2a_input/${name}_h3k4me1_r2.bw  --numberOfProcessors 24
bamCoverage -b $g2a_input/${name_alt}_h3k4me1_r1.bam -o $g2a_input/${name_alt}_h3k4me1_r1.bw  --numberOfProcessors 24


# 4
# DONE
# Get all genes used for the predictions and generate a file with ~22000 lines, one per gene, with 250+/- TSS position.
# This will generate the `regions_includelist` to run ABC's `src/makeCandidateRegions.py`
# chr start end gene 0 strand
awk -v OFS="\t" '{if($6=="-"){$2=$3} {print $1,$2-250,$2+250,$4,$5,$6}}' $g2a/genes_tss_tts.bed > $g2a/whitelist_genes_tss.bed
wc -l $g2a/whitelist_genes_tss.bed

# 5
# DONE
# Transform the RNA-seq expression data for B72 into the TPM-normalized two-column tab-delimited gene-expression format
tmp_expression=`mktemp`
tmp_order=`mktemp`
cut -f1,12 -d, $src/data/rnaseq.csv |tr ',' '\t' | sed 1d | sort -k1,1 > $tmp_expression
sort -k4,4 $g2a/ABC_input/genes_tss_tts.bed > $tmp_order
join $tmp_order $tmp_expression -1 4 -2 1 | tr ' ' '\t' | sort -k2,2 -k3,3n | cut -f1,7 > $g2a/ABC_input/${name}_expression.tsv
cut -f1,11 -d, $src/data/rnaseq.csv |tr ',' '\t' | sed 1d | sort -k1,1 > $tmp_expression
sort -k4,4 $g2a/ABC_input/genes_tss_tts.bed > $tmp_order
join $tmp_order $tmp_expression -1 4 -2 1 | tr ' ' '\t' | sort -k2,2 -k3,3n | cut -f1,7 > $g2a/ABC_input/${name_alt}_expression.tsv
wc -l $g2a/ABC_input/*_expression.tsv
head $g2a/ABC_input/*_expression.tsv

# 6 
# DONE
# Get the 4781 ubiquitously expressed genes:
# A Comprehensive Mouse Transcriptomic BodyMap across 17 Tissues by RNA-seq; Bin Li
# Obtainable from https://www.nature.com/articles/s41598-017-04520-z#Sec2 linked as `Dataset 1`
# 41598_2017_4520_MOESM2_ESM.xls in downloads folder
wc -l $g2a/ABC_input/ubiquitous_expressed_genes.csv

# 7
# DONE
# Get chromosome sizes
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes -O $g2a/ABC_input/mm10.chrom.sizes.unfiltered
sort -k1,1 -k2,2n $g2a/ABC_input/mm10.chrom.sizes.unfiltered > $g2a/ABC_input/mm10.chrom.sizes.extra
grep -v -e Un_ -e _random $g2a/ABC_input/mm10.chrom.sizes.extra > $g2a/ABC_input/mm10.chrom.sizes
rm $g2a/mm10.chrom.sizes.unfiltered
awk -v OFS="\t" '{$2=0"\t"$2; print}' $g2a/ABC_input/mm10.chrom.sizes > $g2a/ABC_input/mm10.chrom.sizes.bed
awk -v OFS="\t" '{$2=0"\t"$2; print}' $g2a/ABC_input/mm10.chrom.sizes.extra > $g2a/ABC_input/mm10.chrom.sizes.extra.bed

# 8
# Reformat Hi-C data to fit ABC's expected input
# bedpe
# chr pos1 pos2 chr pos1 pos2 . signal

# /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/01.HiC_CasellasPublicDatasets/02.Analysis/MergedData/Bcell_72h/hic_results/matrix/Bcell_72h/raw/10000
hic_activ72h_coords=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/01.HiC_CasellasPublicDatasets/02.Analysis/MergedData/Bcell_72h/hic_results/matrix/Bcell_72h/raw/10000/Bcell_72h_10000_abs.bed
hic_activ72h_raw=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/01.HiC_CasellasPublicDatasets/02.Analysis/MergedData/Bcell_72h/hic_results/matrix/Bcell_72h/raw/10000/Bcell_72h_10000.matrix
hic_activ72h_iced=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input/hic/activ72h_B/iced/activ72h_B_10000_iced.matrix.tab

hic_resting_coords=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/01.HiC_CasellasPublicDatasets/02.Analysis/MergedData/Bcell_00h/hic_results/matrix/Bcell_00h/raw/10000/Bcell_00h_10000_abs.bed
hic_resting_raw=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/01.HiC_CasellasPublicDatasets/02.Analysis/MergedData/Bcell_00h/hic_results/matrix/Bcell_00h/raw/10000/Bcell_00h_10000.matrix
hic_resting_iced=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input/hic/${name_alt}/iced/${name_alt}_10000_iced.matrix.tab

# Awk code that filter chr_x chr_x interactions
cd $g2a_hic
awk -v outdir=$g2a_hic/raw 'NR==FNR {  # Process the translation file
        id_to_chr[$4] = $1
        id_to_pos[$4] = $1"\t"$2"\t"$3
        next
}
{  # Process the data file
    if (id_to_chr[$1] == id_to_chr[$2]) {
        output_file = id_to_chr[$1]".bedpe"
        print id_to_pos[$1]"\t"id_to_pos[$2]"\t.\t"$3 >> outdir"/"output_file
    }
}' $hic_activ72h_coords $hic_activ72h_raw 

awk -v outdir=$g2a_hic/iced 'NR==FNR {  # Process the translation file
        id_to_chr[$4] = $1
        id_to_pos[$4] = $1"\t"$2"\t"$3
        next
}
{  # Process the data file
    if (id_to_chr[$1] == id_to_chr[$2]) {
        output_file = id_to_chr[$1]".bedpe"
        print id_to_pos[$1]"\t"id_to_pos[$2]"\t.\t"$3 >> outdir"/"output_file
    }
}' $hic_activ72h_coords $hic_activ72h_iced 

for i in {1..19} Y X M; do
  mkdir -p $g2a_hic/iced/chr${i}
  mv $g2a_hic/iced/chr${i}.bedpe $g2a_hic/iced/chr${i}/
  gzip $g2a_hic/iced/chr${i}/chr${i}.bedpe &
done

# >>> Resting
# Awk code wthat filter chr_x chr_x interactions
cd $g2a_hic_alt
awk -v outdir=$g2a_hic_alt/raw 'NR==FNR {  # Process the translation file
        id_to_chr[$4] = $1
        id_to_pos[$4] = $1"\t"$2"\t"$3
        next
}
{  # Process the data file
    if (id_to_chr[$1] == id_to_chr[$2]) {
        output_file = id_to_chr[$1]".bedpe"
        print id_to_pos[$1]"\t"id_to_pos[$2]"\t.\t"$3 >> outdir"/"output_file
    }
}' $hic_resting_coords $hic_resting_raw 

awk -v outdir=$g2a_hic_alt/iced 'NR==FNR {  # Process the translation file
        id_to_chr[$4] = $1
        id_to_pos[$4] = $1"\t"$2"\t"$3
        next
}
{  # Process the data file
    if (id_to_chr[$1] == id_to_chr[$2]) {
        output_file = id_to_chr[$1]".bedpe"
        print id_to_pos[$1]"\t"id_to_pos[$2]"\t.\t"$3 >> outdir"/"output_file
    }
}' $hic_resting_coords $hic_resting_iced 

for i in {1..19} Y X M; do
#   mkdir -p $g2a_hic/iced/chr${i}
#   mv $g2a_hic/iced/chr${i}.bedpe $g2a_hic/iced/chr${i}/
#   gzip $g2a_hic/iced/chr${i}/chr${i}.bedpe &
  mkdir -p $g2a_hic_alt/iced/chr${i}
  mv $g2a_hic_alt/iced/chr${i}.bedpe $g2a_hic_alt/iced/chr${i}/
  gzip $g2a_hic_alt/iced/chr${i}/chr${i}.bedpe &
done

# mv example hic files to a lotation for observation:
hic_example=`mktemp -d`
abc_example=/mnt/bioadhoc-temp/Groups/RaoLab/Edahi/ForFerhatGit/ABC/example_chr22/input_data/HiC/raw/chr22
cp $abc_example/* $hic_example/
gunzip $hic_example/*
head $hic_example/*
# python src/juicebox_dump.py \
# --hic_file https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined_30.hic \
# --juicebox "java -jar juicer_tools.jar" \
# --outdir /scratch/tmp.6hyFanB8vF/ \
# --chromosomes 22
