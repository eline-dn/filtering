"""
steps:
- load all df in pandas
- add a "conditions" column for each
- apply the filtering function + writes out successful binders + saves modified df
- vertically "concat" them 
- plot iptm pair plots and condition box plots

"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import shutil
import os
from filters_utils import *

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
  design_folder=f"{BC_output}/n{i}/Accepted"
  condition=f"n{i}"
  # load the metrics.csv df
  df=pd.read_csv(csv)
  df['condition']=condition
  print(f"Filtering df for condition {condition}:----------")
  #filtering
  df=compute_specific_metrics(df)
  df=compute_nonspecific_metrics(df)
  ## filters on non_specific reprediction scores
  filtered_df_non_specific = df[(df['mean_empty_i_pTM'] <0.5) & (df['mean_plugged_i_pTM'] >=0.8) & (df['interface_dSASA']>=1700)]
  ## filters on specific reprediction scores
  filtered_df_specific = df[(df['specific_mean_empty_i_pTM'] <0.5) & (df['specific_mean_plugged_i_pTM'] >=0.8)& (df['interface_dSASA']>=1700)]
  ## exctract the interesting binders
  extract_filtered_binders(filtered_df_non_specific, design_folder, filtered_binders_pdbs_non_spec)
  extract_filtered_binders(filtered_df_specific, design_folder, filtered_binders_pdbs_spec)
  
  # concat:
   print("Concatenating to global df--------------")
   global_df=pd.concat([global_df, df])


# boxplots
columns_to_plot = ['interface_dSASA', 'Average_n_InterfaceResidues', 'interface_interface_hbonds', 'interface_hbond_percentage']
condition_boxplots(global_df, columns_to_plot, output_folder)

# ipTM pairplots
columns = ['specific_mean_empty_i_pTM', 'specific_mean_plugged_i_pTM']
_ = sns.pairplot(
    data=global_df,
    vars=columns,
    hue=condition,
    plot_kws={"alpha": 0.2},
    height=3,
    diag_kind="hist",
    diag_kws={"bins": 30},)

plt.savefig(f"{output_folder}/pairplot_specific_iptms.png")

columns = ['mean_empty_i_pTM', 'mean_plugged_i_pTM']
 _ = sns.pairplot(
    data=global_df,
    vars=columns,
    hue=condition,
    plot_kws={"alpha": 0.2},
    height=3,
    diag_kind="hist",
    diag_kws={"bins": 30},)
plt.savefig(f"{output_folder}/pairplot_iptms2.png")



# save global df
global_df.to_csv(f"{output_folder}/global_metrics_df.csv")
