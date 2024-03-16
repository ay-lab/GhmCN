# Sconda conda activate ghmc1; cd /mnt/BioAdHoc/Groups/RaoLab/Edahi/Projects/09.TimeCourse5hmC/ForFerhatGit/GhmCN/reviewer_questions/plot_showing_significance_of_top_10_nearest_10; python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel

# Re-define the data sets as the environment has been reset
font_name = 'Nimbus Roman'

# cells
cells = ['Resting B', 'Activ72h B', 'DP', 'Naive CD4T', 'Naive CD8T', 'Th2']

# Data AUC Score
matched_model_auc = [0.8604, 0.8586, 0.811, 0.8618, 0.855, 0.8198]
nearest_10_auc = [0.8304, 0.8304, 0.7693, 0.8162, 0.8103, 0.7728]
average_model_auc = [0.8134, 0.7941, 0.7875, 0.8423, 0.8599, 0.7897]

# Data AUPR Score
matched_model_aupr = [0.8158, 0.7981, 0.781, 0.8367, 0.8343, 0.7874]
nearest_10_aupr = [0.8124, 0.7878, 0.7488, 0.7778, 0.7773, 0.7603]
average_model_aupr = [0.7815, 0.7649, 0.7552, 0.8113, 0.83, 0.7568]

# Next 5
t_stat, p_value_auc = ttest_rel(matched_model_auc, nearest_10_auc)
t_stat_new, p_value_aupr = ttest_rel(matched_model_aupr, nearest_10_aupr)
# Averaged
t_stat_ave, p_value_auc_ave = ttest_rel(matched_model_auc, average_model_auc)
t_stat_new_ave, p_value_aupr_ave = ttest_rel(matched_model_aupr, average_model_aupr)

# Combining the data into a single DataFrame
df = pd.DataFrame({
    'Top 10 (3D) AUC': matched_model_auc,
    'Nearest 10 (1D) AUC': nearest_10_auc,
    'Top 10 (3D) AUPR': matched_model_aupr,
    'Nearest 10 (1D) AUPR': nearest_10_aupr,
})

df_mix = pd.DataFrame({
    'Matched AUC': matched_model_auc,
    'Average AUC': average_model_auc,
    'Matched AUPR': matched_model_aupr,
    'Average AUPR': average_model_aupr,
})

# Melting the DataFrame
df_melted = df.melt(var_name='Group', value_name='Score')
df_melted_mix = df_mix.melt(var_name='Group', value_name='Score')

# Adding a Model identifier (will deprecate "Pair" if successful)
n_val_per_cell = df_melted.shape[0]//len(cells)
df_melted['Model'] = cells * n_val_per_cell
df_melted_mix['Model'] = cells * n_val_per_cell

# Create a new 'Metric' column that indicates whether the score is AUC or AUPR
df_melted['Metric'] = df_melted['Group'].str.extract('(AUC|AUPR)')
df_melted_mix['Metric'] = df_melted_mix['Group'].str.extract('(AUC|AUPR)')

# Create a 'Condition' column that indicates whether the group is 'Top 10 (3D)' or 'Nearest 10 (1D)'
df_melted['Condition'] = df_melted['Group'].str.extract('(Top 10 \(3D\)|Nearest 10 \(1D\))')
df_melted_mix['Condition'] = df_melted_mix['Group'].str.extract('(Matched|Average)')


def add_significance_bar(ax, y_max, group_col, y_offset=0.02):
    x1, x2 = 0, 1  # Positions of boxplots
    ax.plot([x1, x1, x2, x2], [y_max+0.01, y_max+y_offset, y_max+y_offset, y_max+0.01], lw=1.5, c='k')
    ax.text((x1+x2)*0.5, y_max+y_offset, "*", ha='center', va='bottom', color='k')
    return y_max + y_offset*2  # Return position for p-value text


#
#
#
# # vs Next 5
#
#
#
# Set up the FacetGrid
# plt.close()
g = sns.FacetGrid(df_melted, col="Metric", height=5, aspect=7/10, sharey=True)

# Add the boxplots first
g.map(sns.boxplot, "Condition", "Score", "Group", palette="pastel", order=["Top 10 (3D)", "Nearest 10 (1D)"])

# Now add the point plots for paired lines
# We use 'zorder=10' to draw the points and lines on top of the boxplots
g.map(sns.pointplot, "Condition", "Score", "Model", order=["Top 10 (3D)", "Nearest 10 (1D)"], hue_order = cells,
     palette="viridis", markers="", linestyles="-", linewidth=2.5, alpha=1.0, zorder=10)

max_AUC = df_melted.groupby(['Metric','Condition']).max().Score[( 'AUC')].max()
max_AUPR = df_melted.groupby(['Metric','Condition']).max().Score[( 'AUPR')].max()

for ax, (metric_name, metric_data) in zip(g.axes.flatten(), df_melted.groupby("Metric")):
    ax.set_ylabel('Score', fontsize=18, fontname=font_name)
    ax.set_xlabel('Model', fontsize=16, fontname=font_name)
    for label in ax.get_xticklabels():
        label.set_fontsize(14)
        label.set_fontname(font_name)
    for label in ax.get_yticklabels():
        label.set_fontsize(14)
        label.set_fontname(font_name)
    p_value = p_value_auc if metric_name == "AUC" else p_value_aupr
    y_max = max_AUC if metric_name == "AUC" else max_AUPR
    if p_value < 0.05:
        y_text = add_significance_bar(ax, y_max, 'Condition')
        ax.set_title(f'{metric_name} Score\nPaired t - pval: {p_value:.1e}', fontsize=14, fontname=font_name)

g.fig.subplots_adjust(top=0.9) # adjust the Figure in seaborn

g.add_legend(title='Model', prop={'size': 12, 'family': font_name})
plt.savefig('./merged_significance_dotted_random_split_test_Matched_vs_Next5.pdf')
plt.close()


#
#
#
# # vs Averaged HiC
#
#
#
# Set up the FacetGrid
# plt.close()
g = sns.FacetGrid(df_melted_mix, col="Metric", height=5, aspect=7/10, sharey=True)

# Add the boxplots first
g.map(sns.boxplot, "Condition", "Score", "Group", palette="pastel", order=["Matched", "Average"])

# Now add the point plots for paired lines
# We use 'zorder=10' to draw the points and lines on top of the boxplots
g.map(sns.pointplot, "Condition", "Score", "Model", order=["Matched", "Average"], hue_order = cells,
     palette="viridis", markers="", linestyles="-", linewidth=2.5, alpha=1.0, zorder=10)

max_AUC = df_melted_mix.groupby(['Metric','Condition']).max().Score[( 'AUC')].max()
max_AUPR = df_melted_mix.groupby(['Metric','Condition']).max().Score[( 'AUPR')].max()

for ax, (metric_name, metric_data) in zip(g.axes.flatten(), df_melted_mix.groupby("Metric")):
    ax.set_ylabel('Score', fontsize=18, fontname=font_name)
    ax.set_xlabel('Model', fontsize=16, fontname=font_name)
    for label in ax.get_xticklabels():
        label.set_fontsize(14)
        label.set_fontname(font_name)
    for label in ax.get_yticklabels():
        label.set_fontsize(14)
        label.set_fontname(font_name)
    p_value = p_value_auc_ave if metric_name == "AUC" else p_value_aupr_ave
    y_max = max_AUC if metric_name == "AUC" else max_AUPR
    if p_value < 0.05:
        y_text = add_significance_bar(ax, y_max, 'Condition')
        ax.set_title(f'{metric_name} Score\nPaired t - pval: {p_value:.1e}', fontsize=14, fontname=font_name)

g.fig.subplots_adjust(top=0.9) # adjust the Figure in seaborn

g.add_legend(title='Model', prop={'size': 12, 'family': font_name})

plt.savefig('./merged_significance_dotted_random_split_test_Matched_vs_Average.pdf')
plt.close()
