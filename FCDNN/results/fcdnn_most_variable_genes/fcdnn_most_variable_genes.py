# Sconda conda activate abc; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_most_variable_genes; python
# NOTE: Run /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/results/fcdnn_most_variable_genes/fcdnn_most_variable_genes.sh before this code
import pandas as pd
import numpy as np
import os
# Nupl2   High    1       chr5    1       22.7810339809189

out_dir='/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/'

concatenated_y = None
concatenated_x = None
for subset in ['Dev','Test','Train']:
    m_set = pd.read_csv(out_dir + f'CleanSamples_{subset}_M.txt',names= ['genes','label','binary_label','chr','rank','tpm'],usecols=['genes'],sep="\t")
    m_set.reset_index(drop=False,inplace=True)
    m_set_genes = pd.read_csv(f'/tmp/{subset}.txt',names= ['genes'],sep="\t")
    m_set = m_set.merge(m_set_genes, how='inner',on=['genes'])
    x_data_sub = np.load(f'/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_{subset}_X.npz')['arr_0']
    y_data_sub = np.load(f'/mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/FCDNN/data/C2B_Merged/CleanSamples_{subset}_Y.npz')['arr_0']
    print(f'X orig {subset}: {x_data_sub.shape}')
    print(f'Y orig {subset}: {y_data_sub.shape}')
    assert m_set.shape[0] == x_data_sub[:,m_set.index].shape[1]
    assert m_set.shape[0] == y_data_sub[:,m_set.index].shape[1]
    x_data_sub = x_data_sub[:,m_set.index]
    y_data_sub = y_data_sub[:,m_set.index]
    print(f'X new {subset}: {x_data_sub.shape}')
    print(f'Y new {subset}: {y_data_sub.shape}')
    if concatenated_x is None:
        concatenated_x = x_data_sub
        concatenated_y = y_data_sub
    else:
        concatenated_x = np.concatenate((concatenated_x, x_data_sub), axis=1)
        concatenated_y = np.concatenate((concatenated_y, y_data_sub), axis=1)

np.savez(os.path.join(out_dir, 'CleanSamplesMostVariable_X.npz'), concatenated_x)
np.savez(os.path.join(out_dir, 'CleanSamplesMostVariable_Y.npz'), concatenated_y)

# Test loading:
x_data = np.load(os.path.join(out_dir, 'CleanSamplesMostVariable_X.npz'))['arr_0']
y_data = np.load(os.path.join(out_dir, 'CleanSamplesMostVariable_Y.npz'))['arr_0']
x_data.shape
y_data.shape

