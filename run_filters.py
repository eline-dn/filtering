"""
steps:
- load all df in pandas
- add a "conditions" column for each
- apply the filtering function + writes out successful binders + saves modified df
- vertically "concat" them 
- plot iptm pair plots and condition box plots

"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import shutil
import os
from filters_utils import *
import time

script_start_time = time.time()
global_df=pd.DataFrame()
BC_output=sys.argv[1] # path to the folder where every condition's BC output folder can be found. Each condition's folder must contain an Accepted folder with the binders' pdbs
metrics_output=sys.argv[2] # path to the folder where every condition's run_metrics.py output folder can be found. Each folder has a ./metrics.csv 
output_folder=sys.argv[3] # where to store the successful binders and plots

 
# dealing with 1st reprediction metrics (specific target)
filtered_binders_pdbs_spec=f"{output_folder}/filtered_binders_spec"
os.makedirs(filtered_binders_pdbs_spec, exist_ok=True)
#dealing with metrics from 2nd prediction (aspecific target)
filtered_binders_pdbs_non_spec=f"{output_folder}/filtered_binders_non_spec"
os.makedirs(filtered_binders_pdbs_non_spec, exist_ok=True)

cond_list=[3,9,10,11,12]

for i in cond_list:
  csv=f"{metrics_output}/n{i}/metrics.csv"
  accepted_folder=f"{BC_output}/n{i}/Accepted"
  condition=f"n{i}"
  # load the metrics.csv df
  df=pd.read_csv(csv)
  df['condition']=condition
  print(f"Filtering df for condition {condition}:----------")
  #filtering
  df=compute_specific_metrics(df)
  df=compute_nonspecific_metrics(df)
  df=compute_cofolding_metrics(df)
  ## filters on non_specific reprediction scores
  filtered_df_non_specific = df[(df['mean_empty_i_pTM'] <0.5) & (df['mean_plugged_i_pTM'] >=0.8) & (df['interface_dSASA']>=1700)]
  ## filters on specific reprediction scores
  filtered_df_specific = df[(df['specific_mean_empty_i_pTM'] <0.5) & (df['specific_mean_plugged_i_pTM'] >=0.8)& (df['interface_dSASA']>=1700)]
  ## exctract the interesting binders
  extract_filtered_binders(filtered_df_non_specific, accepted_folder, filtered_binders_pdbs_non_spec)
  extract_filtered_binders(filtered_df_specific, accepted_folder, filtered_binders_pdbs_spec)
  # concat:
  print("Concatenating to global df--------------")
  global_df=pd.concat([global_df, df])

# save global df
global_df.to_csv(f"{output_folder}/global_metrics_df.csv")

# boxplots
columns_to_plot = ['interface_dSASA', 'Average_n_InterfaceResidues', 'interface_interface_hbonds', 'interface_hbond_percentage']
condition_boxplots(global_df, columns_to_plot, output_folder)

# ipTM pairplots
columns = ['specific_mean_empty_i_pTM', 'specific_mean_plugged_i_pTM']
_ = sns.pairplot(
    data=global_df,
    vars=columns,
    hue='condition',
    plot_kws={"alpha": 0.2},
    height=3,
    diag_kind="hist",
    diag_kws={"bins": 30},)

plt.savefig(f"{output_folder}/pairplot_specific_iptms.png")

columns = ['mean_empty_i_pTM', 'mean_plugged_i_pTM']
_ = sns.pairplot(
    data=global_df,
    vars=columns,
    hue='condition',
    plot_kws={"alpha": 0.2},
    height=3,
    diag_kind="hist",
    diag_kws={"bins": 30},)
plt.savefig(f"{output_folder}/pairplot_iptms2.png")



# plot the ipTM scatter plots with thresholds and compare them, highlight selected binders from the first round of selection
selected1_binder_list=global_df[(global_df['specific_mean_empty_i_pTM'] <0.5) & (global_df['specific_mean_plugged_i_pTM'] >=0.8)& (global_df['interface_dSASA']>=1700)].Design
global_df['selected1'] = global_df['Design'].isin(selected1_binder_list)

# scatter plot for case n°1
_=sns.scatterplot(data=global_df, x='specific_mean_plugged_i_pTM', y='specific_mean_empty_i_pTM',hue='selected1')
plugged_ipTM = 0.8
plt.axvline(x=plugged_ipTM, ymin=0, ymax=1, color="black", linestyle="--")
empty_ipTM = 0.5
plt.axhline(
    y=empty_ipTM, xmin=0, xmax=1, color="black", linestyle="--"
)
plt.xlabel("Plugged i_pTM ")
plt.ylabel("Empty i_pTM")
plt.title("ipTM Scatterplot for the ipTM values from the reprediction with target from BC (case 1)")
plt.legend(title="Binders selected with case 1 filters")
plt.savefig("{output_folder}/scatter_1.png")

# scatter plot for case n°2
_=sns.scatterplot(data=global_df, x='mean_plugged_i_pTM', y='mean_empty_i_pTM',hue='selected1')
plugged_ipTM = 0.8
plt.axvline(x=plugged_ipTM, ymin=0, ymax=1, color="black", linestyle="--")
empty_ipTM = 0.5
plt.axhline(
    y=empty_ipTM, xmin=0, xmax=1, color="black", linestyle="--"
)
plt.xlabel("Plugged i_pTM ")
plt.ylabel("Empty i_pTM")
plt.title("ipTM Scatterplot for the ipTM values from the reprediction with initial target (case 2)")
plt.legend(title="Binders selected with case 1 filters")
plt.savefig("{output_folder}/scatter_2.png")

# scatter plot for case n°3
_=sns.scatterplot(data=global_df, x='cofold_mean_plugged_i_pTM', y='cofold_mean_empty_i_pTM',hue='selected1')
plugged_ipTM = 0.8
plt.axvline(x=plugged_ipTM, ymin=0, ymax=1, color="black", linestyle="--")
empty_ipTM = 0.5
plt.axhline(
    y=empty_ipTM, xmin=0, xmax=1, color="black", linestyle="--"
)
plt.xlabel("Plugged i_pTM ")
plt.ylabel("Empty i_pTM")
plt.title("ipTM Scatterplot for the ipTM values from the reprediction without target template (case 3)")
plt.legend(title="Binders selected with case 1 filters")
plt.savefig("{output_folder}/scatter_3.png")



# end of the script: how long?
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
print("Finished filtering. Script execution took: "+elapsed_text)
