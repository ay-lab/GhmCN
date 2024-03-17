# Sconda conda activate plot_gc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_as_expression_predictor; python
import pybedtools
import pyfaidx
import pandas as pd

def convert_tpm_to_labels(df, balanced_split=False):
    # Calculate mean per column, excluding zeroes
    medians = df.replace(0, pd.NA).median() if not balanced_split else df.median()
    # Update values to 'above' or 'below'
    for column in df.columns:
        df[column] = df[column].apply(lambda x: 1 if x > medians[column] else 0)
    return df

def calculate_cpg_content(sequence):
    cpg_count = sequence.count('CG')
    total_count = len(sequence)/2
    cpg_content = (cpg_count / total_count) * 100
    return cpg_content

def get_cpg_content(bed_file,genome_fasta):
    # Load genome fasta file
    genome = pyfaidx.Fasta(genome_fasta)
    peaks = pybedtools.BedTool(bed_file)
    cpg_contents = []
    # Calculate GC content for each peak
    for peak in peaks:
        sequence = str(genome[peak.chrom][peak.start:peak.end]).upper()
        if 'N' in  sequence:
            cpg_contents.append(0)
            continue
        cpg_content = calculate_cpg_content(sequence)
        cpg_contents.append(cpg_content)
    return cpg_contents



out_path='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_as_expression_predictor'
bed_header_format = ['chr', 'start', 'end', 'gene']
# Load all promoters
promoter_coords = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences/all_genes_promoter_coords.bed',sep="\t")
# set genome
genome_fasta = '/mnt/BioAdHoc/Groups/RaoLab/Edahi/01.Genomes/mm10/mm10.fa'

bed_file = 'promoter_coords_of_genes_with_label.bed'
cpg_contents = get_cpg_content(bed_file, genome_fasta)



# # New split
# Load all sample's calculated TPM
all_samples_tpm = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences/FILTERED_all_samples_tpm.csv')

# Update TPMs df to labels
all_samples_tpm.iloc[:,1::]  =  convert_tpm_to_labels(all_samples_tpm.drop(columns = {'gene'}))

# Subset to genes with promoters
all_samples_tpm = all_samples_tpm.merge(promoter_coords, on = ['gene'], how = 'inner')
all_samples_tpm[bed_header_format].to_csv('FILTERED_promoter_coords_of_genes_with_label.bed', index=False, header=False, sep="\t")

all_samples_tpm['cpg_content'] = cpg_contents
all_samples_tpm.drop(columns = ['chr','start','end']).to_csv('FILTERED_all_samples_classes_and_cpg_content_in_promoter_new_split.tsv', index=False, sep="\t")

# # Balanced split
# Load all sample's calculated TPM
all_samples_tpm = pd.read_csv('/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/cpg_composition_differences/FILTERED_all_samples_tpm.csv')

# Update TPMs df to labels
all_samples_tpm.iloc[:,1::]  =  convert_tpm_to_labels(all_samples_tpm.drop(columns = {'gene'}), balanced_split=True)

# Subset to genes with promoters
all_samples_tpm = all_samples_tpm.merge(promoter_coords, on = ['gene'], how = 'inner')

all_samples_tpm['cpg_content'] = cpg_contents
all_samples_tpm.drop(columns = ['chr','start','end']).to_csv('FILTERED_all_samples_classes_and_cpg_content_in_promoter_balanced_split.tsv', index=False, sep="\t")
