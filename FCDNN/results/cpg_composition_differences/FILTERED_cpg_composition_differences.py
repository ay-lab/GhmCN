# Sconda conda activate plot_gc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences; python
import pybedtools
import pyfaidx
import matplotlib.pyplot as plt
import pandas as pd
import os


# Function to TPM normalize the gene counts
def tpm_normalize(df):
    # Assuming column 6 has gene length in base pairs (bp)
    length_kb = df.iloc[:, 5] / 1000  # Convert gene length from bp to kilobases (kb)
    # Assuming column 7 has raw read counts
    reads = df.iloc[:, 6]
    # Calculate Reads Per Kilobase (RPK)
    rpk = reads / length_kb
    # Calculate the scaling factor by summing all RPK values and dividing by 1,000,000
    scaling_factor = rpk.sum() / 1_000_000
    # Divide RPK by the scaling factor to get TPM
    tpm = rpk / scaling_factor
    return tpm


# Add promoter +/- 1500
def get_promoter_start(row, extend:int = 1500):
    promoter_idx = 'start' if row['direction'] == "+" else 'end'
    start = row[promoter_idx] - extend
    return start


def get_promoter_end(row, extend:int = 1500):
    promoter_idx = 'start' if row['direction'] == "+" else 'end'
    end = row[promoter_idx] + extend
    return end


def calculate_cpg_content(sequence):
    cpg_count = sequence.count('CG')
    total_count = len(sequence)/2
    cpg_content = (cpg_count / total_count) * 100
    return cpg_content


# Function to plot GC content for each peak set
def plot_cpg_content(bed_files:list, genome_fasta:str, titles:list, transparency:list=None, colors=('skyblue', 'red'),log=False, density=False):
    # Load genome fasta file
    genome = pyfaidx.Fasta(genome_fasta)
    logname=""
    densityname=""
    transparency = [0.7]*len(bed_files) if transparency is None else transparency
    all_cpg_contents = []
    # Plotting with color and transparency
    for i, bed_file in enumerate(bed_files):
        peaks = pybedtools.BedTool(bed_file)
        cpg_contents = []
        
        # Calculate GC content for each peak
        for peak in peaks:
            sequence = str(genome[peak.chrom][peak.start:peak.end]).upper()
            if 'N' in  sequence:
                continue
            cpg_content = calculate_cpg_content(sequence)
            cpg_contents.append(cpg_content)
        all_cpg_contents.append(cpg_contents)
        # Plotting for each bed file
        plt.hist(cpg_contents,density=density, bins=20, color=colors[i], edgecolor='black', alpha=transparency[i], label=titles[i])
    
    promoter_size = peaks.to_dataframe().loc[0,:]
    promoter_size = abs(promoter_size['end'] - promoter_size['start'])
    plt.xlabel('CpG Content (%)')
    ylabel='Frequency'
    if log:
        plt.yscale('log')  # Set y-axis to log scale
        logname="_log"
        ylabel+=" (log10)"
    if density:
        densityname = "_density"
        ylabel+=" (density)"
    plt.ylabel(ylabel)
    plt.title(f'CpG Content Distribution\n{promoter_size}bp Promoter')
    plt.legend()
    
    # Save the plot as a PDF in the working directory
    plt.savefig(os.path.join(os.getcwd(), f'FILTERED_CpG_Content_Distribution_{promoter_size}bp_promoter_size{logname}{densityname}_noN.pdf'))
    
    # Show the plot
    plt.close()
        
    return(all_cpg_contents)


out_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences'
bed_header_format = ['chr', 'start', 'end', 'gene']
# Read the file with rna files location
rna_files_location = pd.read_csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/filtered_files_surnames.csv", header = None, usecols=[7,10])
rna_files_location.columns = ['sample', 'rna_file']

# Get all data's TPM and save matrix
if not os.path.exists(os.path.join(out_path, "FILTERED_all_samples_tpm.csv")):
    all_samples_tpm = []
    for sample in rna_files_location['sample']:
        rna_file = rna_files_location.query('sample == @sample')["rna_file"].reset_index(drop=True)[0]
        rna_data = pd.read_csv(rna_file, sep='\t', comment='#')
        all_samples_tpm.append(tpm_normalize(rna_data).tolist())
        print(sample)
    
    all_samples_tpm = pd.DataFrame(all_samples_tpm).T
    all_samples_tpm = pd.concat([rna_data.Geneid.copy(), all_samples_tpm], axis=1)
    all_samples_tpm.columns = ['gene'] + rna_files_location['sample'].tolist()
    
    # Save/Load all sample's calculated TPM
    all_samples_tpm.to_csv(os.path.join(out_path, "FILTERED_all_samples_tpm.csv"), index=False)

all_samples_tpm = pd.read_csv(os.path.join(out_path, "FILTERED_all_samples_tpm.csv"))

# Get the gene locations
gene_file = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input/genes_tss_tts.bed'
gene_coords = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input/genes_tss_tts.bed', sep="\t", header=None)
gene_coords.columns = ['chr','start','end','gene','x','direction']

# Apply the functions to each row in the DataFrame
promoter_ext = 500
gene_coords['promoter_start'] = gene_coords.apply(get_promoter_start, axis=1, extend = promoter_ext)
gene_coords['promoter_end'] = gene_coords.apply(get_promoter_end, axis=1, extend = promoter_ext)
promoter_coords = gene_coords[['gene','chr','promoter_start','promoter_end']].rename(columns = {'promoter_start':'start', 'promoter_end':'end'})

# # # # # # # #
genome_fasta = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/mm10/mm10.fa'

# Get bed files for
# 1. Ubiquitous genes
ubiquitous_file = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/ghmc_to_abc/ABC_input/ubiquitous_expressed_genes.csv'
ubiquitous_genes = pd.read_csv(ubiquitous_file, header = None, names = ['gene'])
ubiquitous_genes_coords = ubiquitous_genes.merge(promoter_coords, on = ['gene'], how = 'inner')

# 2. Gene expressed / Non expressed across
max_expression_across_samples = all_samples_tpm.drop('gene', axis=1).max(axis=1)
genes_always_off = pd.DataFrame(all_samples_tpm[max_expression_across_samples == 0.0]['gene'])
genes_always_off_coords = genes_always_off.merge(promoter_coords, on = ['gene'], how = 'inner')

# 3. Genes as higher / lower always
path_in = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B'
colnames = ['gene','category','digit_category','chrom','index','tpm']
samples = pd.read_csv("/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/filtered_files_surnames.csv", header = None)[7]
noisy_genes = ["Mir5098", "Eno1b", "Mir684-1", "Gm5643", "Gm5801", "Bc1", "Gm5512", "Btg3", "Mir3473a"]

# Load Sample's gene labels
if not os.path.exists(os.path.join(out_path, "FILTERED_all_classes_49_samples_balanced.csv")):
    for i, name in enumerate(samples.tolist()):
        dev_m = pd.read_table(f'{path_in}/{name}_Dev_M.txt', names=colnames)[['gene', 'category']]
        tst_m = pd.read_table(f'{path_in}/{name}_Test_M.txt', names=colnames)[['gene', 'category']]
        trn_m = pd.read_table(f'{path_in}/{name}_Train_M.txt', names=colnames)[['gene', 'category']]
        all = pd.concat([dev_m, tst_m, trn_m]).rename(columns = {'category':name})
        all = all[~all['gene'].isin(noisy_genes)]
        all_classes = all if i == 0 else pd.merge(all_classes, all, on=["gene"], how="left")
        print(i, name)
    all_classes.to_csv(os.path.join(out_path, "FILTERED_all_classes_49_samples_balanced.csv"), index=False)
else:
    all_classes = pd.read_csv(os.path.join(out_path, "FILTERED_all_classes_49_samples_balanced.csv"))

# Get genes with always same label:
genes_getting_same_label = all_classes[all_classes.drop('gene',axis=1).nunique(axis=1) == 1].iloc[:,[0,1]].iloc[:,[0,1]].reset_index(drop=True)
genes_getting_same_label.columns = ['gene','label']
genes_higher = pd.DataFrame(genes_getting_same_label.query('label == "High"').gene)
genes_lower = pd.DataFrame(genes_getting_same_label.query('label == "Low"').gene)
genes_always_higher_coords = genes_higher.merge(promoter_coords, on = ['gene'], how = 'inner')
genes_always_lower_coords = genes_lower.merge(promoter_coords, on = ['gene'], how = 'inner')

# Get genes with highly variable labels:
# represented as genes whose "High" "Low" distribution across all samples 
# is such that the "Highs" are similar in amount to the "Low"s per gene.
# Meaning that if labels were represented as 1 and -1, then 
# the sum across columns would be closer to 0.
# Since we have 49 samples, an abs number of <15 was considered enough.
genes_any_variable = pd.DataFrame(all_classes.gene[all_classes.drop('gene',axis=1).replace({'High':1, 'Low':-1}).sum(axis=1).abs() < 49]).reset_index(drop=True)
genes_any_variable = genes_any_variable.merge(promoter_coords, on = ['gene'], how = 'inner')
genes_variable = pd.DataFrame(all_classes.gene[all_classes.drop('gene',axis=1).replace({'High':1, 'Low':-1}).sum(axis=1).abs() < 15]).reset_index(drop=True)
genes_variable = genes_variable.merge(promoter_coords, on = ['gene'], how = 'inner')

# Save all BED files
bed_header_format = ['chr', 'start', 'end', 'gene']
prom_size=promoter_ext*2
promoter_coords[bed_header_format].to_csv(f'all_genes_promoter_coords_{prom_size}bp.bed', index=False, sep="\t")
ubiquitous_genes_coords[bed_header_format].to_csv(f'ubiquitous_genes_promoter_coords_{prom_size}bp.bed', index=False, sep="\t",header=False)
genes_always_off_coords[bed_header_format].to_csv(f'always_off_genes_promoter_coords_{prom_size}bp.bed', index=False, sep="\t",header=False)
genes_always_higher_coords[bed_header_format].to_csv(f'FILTERED_always_higher_genes_promoter_coords_{prom_size}bp.bed', index=False, sep="\t", header=False)
genes_always_lower_coords[bed_header_format].to_csv(f'FILTERED_always_lower_genes_promoter_coords_{prom_size}bp.bed', index=False, sep="\t", header=False)
genes_variable[bed_header_format].to_csv(f'FILTERED_variable_label_genes_promoter_coords_{prom_size}bp.bed', index=False, sep="\t", header=False)
genes_any_variable[bed_header_format].to_csv(f'FILTERED_any_degree_of_variable_label_genes_promoter_coords_{prom_size}bp.bed', index=False, sep="\t", header=False)

ubiquitous_genes_coords.shape     # 4378 
genes_always_higher_coords.shape  # 5459 
genes_variable.shape              # 2447 
genes_always_lower_coords.shape   # 3189 
genes_always_off_coords.shape     # 643  

# Example usage with two BED files and custom colors
for prom_size in [1000, 2000, 3000]:
    beds = [
        f'ubiquitous_genes_promoter_coords_{prom_size}bp.bed',
        f'FILTERED_always_higher_genes_promoter_coords_{prom_size}bp.bed',
        f'FILTERED_variable_label_genes_promoter_coords_{prom_size}bp.bed',
        f'FILTERED_always_lower_genes_promoter_coords_{prom_size}bp.bed',
        f'always_off_genes_promoter_coords_{prom_size}bp.bed',
    ]
    for log in [True, False]:
        for density in [True, False]:
            cpg_contents = plot_cpg_content(beds,
                            genome_fasta,
                            titles=[
                                'Ubiquitous Genes (4378)',
                                'Always Higher Label (5459)',
                                '[Most] Variable Label (2447)',
                                'Always Lower Label (3189)',
                                'TPM = 0 (643)',
                                ],
                            colors=('red', 'gray','yellow', 'green','blue'),
                            transparency=[0.7, 0.65, 0.6, 0.55, 0.5],
                            log=log,
                            density=density)


