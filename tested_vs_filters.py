import sys
import numpy as np
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
#BC_output=sys.argv[1] # path to the folder where every condition's BC output folder can be found. Each condition's folder must contain an Accepted folder with the binders' pdbs
metrics_output=sys.argv[1] # path to the ./metrics.csv folder
output_folder=sys.argv[2] # where to store plots?


csv=f"{metrics_output}/metrics.csv"
  # load the metrics.csv df
df=pd.read_csv(csv)

  #filtering
def compute_nonspecific_metrics(df):
  """
  compute models mean for empty and plugged ipTM from the reprediction with non-specific target
  """
  df['mean_it_i_pTM'] = df[['1_it_i_pTM', '2_it_i_pTM']].mean(axis=1)
  df['mean_it_Binder_RMSD_to_binding_site'] = df[['1_it_Binder_RMSD_to_binding_site', '2_it_Binder_RMSD_to_binding_site']].mean(axis=1)

  return df
df=compute_nonspecific_metrics(df)

  ## filters on non_specific reprediction scores AND on specific reprediction scores (ipTM from BC)
#filtered_df = df[(df['mean_it_i_pTM'] >=0.8) & (df['Average_dSASA']>=1700) & (df['Average_i_pTM'] >=0.8)]

# plot final selected binders 
selected12_binder_list=df[(df['mean_it_i_pTM'] >=0.8) & (df['Average_dSASA']>=1700) & (df['Average_i_pTM'] >=0.8)].DesignName
df['selected12'] = df['DesignName'].isin(selected12_binder_list)
# drop binders without binding information:
df = df.dropna(subset=["Binding"])
# convert to bool
df["Binding"] = df["Binding"].astype(bool)
# define category:
conditions = [
    (df["Binding"] == True) & (df["selected12"] == False),
    (df["Binding"] == False) & (df["selected12"] == True),
    (df["Binding"] == True) & (df["selected12"] == True),
    (df["Binding"] == False) & (df["selected12"] == False),
]
choices = ["binding", "pass_filters", "both", "else"]

# Create new column
df["category"] = np.select(conditions, choices, default="unknown")

# scatter plot: binding binders iptM
_=sns.scatterplot(data=df, x='mean_it_i_pTM', y='Average_i_pTM',hue='category')
initial_target_ipTM = 0.8
plt.axvline(x=initial_target_ipTM, ymin=0, ymax=1, color="black", linestyle="--")
BC_target_ipTM = 0.8
plt.axhline(
    y=BC_target_ipTM, xmin=0, xmax=1, color="black", linestyle="--"
)
plt.xlabel("Initial Target i_pTM (case 2)")
plt.ylabel("BC target i_pTM (case 1")
plt.title("Binding Binders vs Filters")
plt.legend(title="Binders selected with case 1 + case 2 filters")
plt.savefig(f"{output_folder}/scatter_bindersvsfilters.png")
plt.close()


df.to_csv(f"{output_folder}/cleaned_df.csv")

# end of the script: how long?
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
print("Finished filtering. Script execution took: "+elapsed_text)
