# overlap of edges between aicda and GhmC
# srun --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=50g --time=04:00:00 --pty bash -i

name=activ72h_B
g2a=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc
g2a_input=$g2a/ABC_input
char_out=$g2a/post_analysis/characterization_abc_k27ac_5hmc
abc_bed_k27=$g2a/ABC_output/Predictions/EnhancerPredictions.bed
abc_bed_hmc=$g2a/ABC_output_5hmc/Predictions/EnhancerPredictions.bed
abc_bed_hmc_resting=$g2a/ABC_output_5hmc_resting_B/Predictions/EnhancerPredictions.bed

k27_unique=$char_out/regions_unique_for_k27.bed
k27_shared=$char_out/regions_shared_for_k27.bed
hmc_unique=$char_out/regions_unique_for_hmc.bed
hmc_shared=$char_out/regions_shared_for_hmc.bed

Sconda 
conda activate abc

bedtools intersect -a ${abc_bed_k27} -b ${abc_bed_hmc} -v | sort -k1,1V -k2,2n | uniq  > $k27_unique
bedtools intersect -a ${abc_bed_k27} -b ${abc_bed_hmc} -u | sort -k1,1V -k2,2n | uniq  > $k27_shared

bedtools intersect -b ${abc_bed_k27} -a ${abc_bed_hmc} -v | sort -k1,1V -k2,2n | uniq  > $hmc_unique
bedtools intersect -b ${abc_bed_k27} -a ${abc_bed_hmc} -u | sort -k1,1V -k2,2n | uniq  > $hmc_shared
wc -l  $k27_unique $k27_shared $abc_bed_k27 $hmc_unique $hmc_shared $abc_bed_hmc

echo 'chr6 122290000 122300000 long1
chr6 122390000 122400000 long2' | tr ' ' '\t' > /tmp/aicda_enh.bed
cat /tmp/aicda_enh.bed
# narrowed down to the peak summit:
echo 'chr6 122293509 122294342 narrow1
chr6 122393397 122393996 narrow2' | tr ' ' '\t' > /tmp/aicda_enh_narrow.bed

bedtools intersect -a /tmp/aicda_enh.bed -b $k27_unique -wb
bedtools intersect -a /tmp/aicda_enh.bed -b $hmc_unique -wb
bedtools intersect -a /tmp/aicda_enh.bed -b $k27_shared -wb
bedtools intersect -a /tmp/aicda_enh.bed -b $hmc_shared -wb

# check if present in resting ABC regions:
bedtools intersect -a /tmp/aicda_enh.bed -b $abc_bed_hmc_resting -wb
bedtools intersect -a /tmp/aicda_enh_narrow.bed -b $abc_bed_hmc_resting -wb

wc -l $k27_unique \
 $k27_shared \
 $hmc_unique \
 $hmc_shared

#    1448 regions_unique_for_k27.bed
#   11874 regions_shared_for_k27.bed
#   17487 regions_unique_fo_hmc.bed
#   11874 regions_shared_fo_hmc.bed

# 4. Get signal in unique regions.
# Sconda
conda deactivate
conda activate deeptools
labs=(hmC atac h3k27ac)
bws=($g2a_input/${name}_5hmc_r1.bw $g2a_input/${name}_atac_rep1.bw $g2a_input/${name}_h3k27ac.bw)
hmc_u=unique_hmc_1kb
k27_u=unique_k27_1kb
hmc_s=shared_hmc_1kb
k27_s=shared_k27_1kb

# > Matrices
# # Unique
# Unique - 5hmCN
computeMatrix reference-point \
 --referencePoint center \
 --numberOfProcessors 24 \
 -a 500 -b 500 \
 -bs 10 \
 --samplesLabel ${labs[@]} \
 --missingDataAsZero \
 -S ${bws[@]} \
 -R $hmc_unique \
 -out $char_out/matrix_${hmc_u}.tab.gz

# Unique - k27ac
computeMatrix reference-point \
 --referencePoint center \
 --numberOfProcessors 24 \
 -bs 10 \
 --samplesLabel ${labs[@]} \
 --missingDataAsZero \
 -S ${bws[@]} \
 -a 500 -b 500 \
 -R $k27_unique \
 -out $char_out/matrix_${k27_u}.tab.gz

# # Shared
# Shared - 5hmCN
computeMatrix reference-point \
 --referencePoint center \
 --numberOfProcessors 24 \
 -a 500 -b 500 \
 -bs 10 \
 --samplesLabel ${labs[@]} \
 --missingDataAsZero \
 -S ${bws[@]} \
 -R $hmc_shared \
 -out $char_out/matrix_${hmc_s}.tab.gz

# Shared - k27ac
computeMatrix reference-point \
 --referencePoint center \
 --numberOfProcessors 24 \
 -a 500 -b 500 \
 -bs 10 \
 --samplesLabel ${labs[@]} \
 --missingDataAsZero \
 -S ${bws[@]} \
 -R $k27_shared \
 -out $char_out/matrix_${k27_s}.tab.gz


# > Plots
# Unique - 5hmCN
ymax=(30 30 40)
zmax=90
plotHeatmap \
 -m $char_out/matrix_${hmc_u}.tab.gz \
 -out $char_out/${hmc_u}.pdf \
 --outFileSortedRegions $char_out/sorted_plotted_regions_${hmc_u}.bed \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax ${ymax[@]} \
 --zMax $zmax \
 --refPointLabel "" \
 --regionsLabel "Signal" \
 --plotTitle ""

# Unique - k27ac
plotHeatmap \
 -m $char_out/matrix_${k27_u}.tab.gz \
 -out $char_out/${k27_u}.pdf \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax ${ymax[@]} \
 --zMax $zmax \
 --refPointLabel "" \
 --regionsLabel "Signal" \
 --plotTitle ""

# # Shared
# Shared - 5hmCN
plotHeatmap \
 -m $char_out/matrix_${hmc_s}.tab.gz \
 -out $char_out/${hmc_s}.pdf \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax ${ymax[@]} \
 --zMax $zmax \
 --refPointLabel "" \
 --regionsLabel "Signal" \
 --plotTitle ""

# Shared - k27ac
plotHeatmap \
 -m $char_out/matrix_${k27_s}.tab.gz \
 -out $char_out/${k27_s}.pdf \
 --outFileSortedRegions $char_out/sorted_plotted_regions_${k27_s}.bed \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax ${ymax[@]} \
 --zMax $zmax \
 --refPointLabel "" \
 --regionsLabel "Signal" \
 --plotTitle ""

# Increased flanking:
# # Shared
# Shared - 5hmCN
plotHeatmap \
 -m $char_out/matrix_${hmc_s}_plus5k_2.tab.gz \
 -out $char_out/${hmc_s}_plus5k_2.pdf \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax ${ymax[@]} \
 --zMax $zmax \
 --refPointLabel 5hmC \
 --regionsLabel "5hmC shared enhancers" \
 --plotTitle 'Epigenetic signal' 

# Find where the long1 and long2 enhancers are located at in the matrices
# long2 is exclusive for 5hmC
echo 'chr6 122290000 122300000 long1
chr6 122390000 122400000 long2' | tr ' ' '\t' > /tmp/aicda_enh.bed
bedtools intersect -a /tmp/aicda_enh.bed -b $char_out/sorted_plotted_regions_${hmc_u}.bed -wb | cut -f8
# gives chr6:122393471-122393971
grep chr6:122393471-122393971 $char_out/sorted_plotted_regions_${hmc_u}.bed -n
wc -l $char_out/sorted_plotted_regions_${hmc_u}.bed
# It is located in row 5948/17487
# That is about a third way down through the plot: 0.340138389


# long1 is present in the shared plot:
bedtools intersect -a /tmp/aicda_enh.bed -b $char_out/sorted_plotted_regions_${k27_s}.bed -wb | cut -f8
# gives chr6:122293644-122294804
#       chr6:122291244-122292612
#       chr6:122297395-122298117
grep chr6:122293644-122294804 $char_out/sorted_plotted_regions_${k27_s}.bed -n
grep chr6:122291244-122292612 $char_out/sorted_plotted_regions_${k27_s}.bed -n
grep chr6:122297395-122298117 $char_out/sorted_plotted_regions_${k27_s}.bed -n
wc -l $char_out/sorted_plotted_regions_${k27_s}.bed
# It is located in row 1701/11874 *
#                      7571/11874
#                      9564/11874
# That is about a seventh     way down through the plot: 0.1432
# That is about a two thirds  way down through the plot: 0.6376
# That is about a four fifths way down through the plot: 0.8054
# 27 1/8 iinchs


# long1 is present in the shared plot:

# narrowed down to:
echo 'chr6 122293509 122294342 narrow1
chr6 122393397 122393996 narrow2' | tr ' ' '\t' > /tmp/aicda_enh_narrow.bed

# Add some more plots for some other marks:
alt_input=/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/post_analysis/characterization_abc_5hmc_00_72h_B
ls $alt_input/${name}_H3*bw

labs=(H3K4me1 H3K4me2)
bws=(
    $alt_input/${name}_H3K4me1.bw
    $alt_input/${name}_H3K4me2.bw
    )
hmc_u=unique_hmc_1kb
k27_u=unique_k27_1kb
k27_s=shared_k27_1kb

# > Matrices
# Unique - 5hmCN
computeMatrix reference-point \
 --referencePoint center \
 --numberOfProcessors 24 \
 -a 500 -b 500 \
 -bs 10 \
 --samplesLabel ${labs[@]} \
 --missingDataAsZero \
 -S ${bws[@]} \
 -R $hmc_unique \
 -out $char_out/exploratory_matrix_${hmc_u}.tab.gz &

# Unique - k27ac
computeMatrix reference-point \
 --referencePoint center \
 --numberOfProcessors 24 \
 -bs 10 \
 --samplesLabel ${labs[@]} \
 --missingDataAsZero \
 -S ${bws[@]} \
 -a 500 -b 500 \
 -R $k27_unique \
 -out $char_out/exploratory_matrix_${k27_u}.tab.gz &

# Shared - k27ac
computeMatrix reference-point \
 --referencePoint center \
 --numberOfProcessors 24 \
 -a 500 -b 500 \
 -bs 10 \
 --samplesLabel ${labs[@]} \
 --missingDataAsZero \
 -S ${bws[@]} \
 -R $k27_shared \
 -out $char_out/exploratory_matrix_${k27_s}.tab.gz &


# > Plots
# Unique - 5hmCN
ymax=100
zmax=100
plotHeatmap \
 -m $char_out/exploratory_matrix_${hmc_u}.tab.gz \
 -out $char_out/exploratory_${hmc_u}.pdf \
 --outFileSortedRegions $char_out/exploratory_sorted_plotted_regions_${hmc_u}.bed \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax $ymax \
 --zMax $zmax \
 --refPointLabel PeakSummit \
 --regionsLabel "5hmC unique" \
 --plotTitle 'ABC 5hmC -v- H3K27Ac' 

# Unique - k27ac
plotHeatmap \
 -m $char_out/exploratory_matrix_${k27_u}.tab.gz \
 -out $char_out/exploratory_${k27_u}.pdf \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax $ymax \
 --zMax $zmax \
 --refPointLabel PeakSummit \
 --regionsLabel "H3K27Ac unique" \
 --plotTitle 'ABC 5hmC -v- H3K27Ac' 

# Shared - k27ac
plotHeatmap \
 -m $char_out/exploratory_matrix_${k27_s}.tab.gz \
 -out $char_out/exploratory_${k27_s}.pdf \
 --outFileSortedRegions $char_out/exploratory_sorted_plotted_regions_${k27_s}.bed \
 --dpi 400 \
 --heatmapHeight 15  \
 --yMax $ymax \
 --zMax $zmax \
 --refPointLabel PeakSummit \
 --regionsLabel "Shared" \
 --plotTitle 'ABC 5hmC -v- H3K27Ac' 
 