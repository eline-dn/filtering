"""
steps:
- load all df in pandas
- add a "conditions" column for each
- vertically "concat" them 
- apply the filtering and plotting iptms function + writes out successful binders + saves modified df
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import shutil
import os
from filters_utils import *

global_df=pd.DataFrame()

cond_list=[3,9,10,11,12]

for i in cond_list:
  metrics_path[f"n{i}"]=f"/n{i}/metrics.csv"


for condition, csv in metrics_path.items():
  df=pd.read_csv(csv)
  df['condition']=condition
  global_df=pd.concat([global_df, df])

output_folder=

# boxplots
columns_to_plot = ['interface_dSASA', 'Average_n_InterfaceResidues', 'interface_interface_hbonds', 'interface_hbond_percentage']
condition_boxplots(global_df, columns_to_plot, output_folder)




