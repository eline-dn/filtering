import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import shutil
import os

def condition_boxplots(global_df, columns_to_plot, output_folder):
  """
  columns_to_plot: list of cols to plot, 
  e.g. columns_to_plot = ['interface_dSASA', 'Average_n_InterfaceResidues', 'interface_interface_hbonds', 'interface_hbond_percentage']
  output_folder: where to save the plot (png)
  just plot a bunch of boxplots to compare the effect of the conditions on some binder's properties
  """
  fig, axes = plt.subplots(nrows=len(columns_to_plot), ncols=1, figsize=(10, 15))

# Flatten the axes array if there's only one row or column for easy iteration
  if len(columns_to_plot) == 1:
      axes = [axes]

# Plot each column on a separate subplot
  for i, col in enumerate(columns_to_plot):
      sns.boxplot(data=global_df, x='condition', y=col, ax=axes[i])
      axes[i].set_title(f'Boxplot of {col} by Condition')
      axes[i].set_xlabel('Condition')
      axes[i].set_ylabel(col)

# Adjust layout to prevent overlapping titles/labels
  plt.tight_layout()
  plt.show()
  plt.savefig(f"{output_folder}/box_plots_conditions.png")

### some helpers functions
def compute_specific_metrics(df):
  """
  compute models mean for specific empty and plugged ipTM
  
  """
  df['specific_mean_empty_i_pTM'] = df[['specific_1_empty_i_pTM', 'specific_2_empty_i_pTM']].mean(axis=1)
  df['specific_mean_plugged_i_pTM'] = df[['specific_1_plugged_i_pTM', 'specific_2_plugged_i_pTM']].mean(axis=1)
  df['specific_mean_plugged_Binder_RMSD_to_binding_site'] = df[['specific_1_plugged_Binder_RMSD_to_binding_site', 'specific_2_plugged_Binder_RMSD_to_binding_site']].mean(axis=1)
  df['specific_mean_empty_Binder_RMSD_to_binding_site'] = df[['specific_1_empty_Binder_RMSD_to_binding_site', 'specific_2_empty_Binder_RMSD_to_binding_site']].mean(axis=1)
  """
  columns = ['specific_mean_empty_i_pTM', 'specific_mean_plugged_i_pTM']
  _ = sns.pairplot(
    data=df,
    vars=columns,
    #hue=target_column,
    plot_kws={"alpha": 0.2},
    height=3,
    diag_kind="hist",
    diag_kws={"bins": 30},)

  plt.savefig(f"{output_folder}/pairplot_specific_iptms.png")
  """
  return df


def compute_nonspecific_metrics(df):
  """
  compute models mean for empty and plugged ipTM from the reprediction with non-specific target
  """
  df['mean_empty_i_pTM'] = df[['1_empty_i_pTM', '2_empty_i_pTM']].mean(axis=1)
  df['mean_plugged_i_pTM'] = df[['1_plugged_i_pTM', '2_plugged_i_pTM']].mean(axis=1)
  df['mean_plugged_Binder_RMSD_to_binding_site'] = df[['1_plugged_Binder_RMSD_to_binding_site', '2_plugged_Binder_RMSD_to_binding_site']].mean(axis=1)
  df['mean_empty_Binder_RMSD_to_binding_site'] = df[['1_empty_Binder_RMSD_to_binding_site', '2_empty_Binder_RMSD_to_binding_site']].mean(axis=1)

  """ 
  columns = ['mean_empty_i_pTM', 'mean_plugged_i_pTM']
  _ = sns.pairplot(
    data=df,
    vars=columns,
    #hue=target_column,
    plot_kws={"alpha": 0.2},
    height=3,
    diag_kind="hist",
    diag_kws={"bins": 30},)
  plt.savefig(f"{output_folder}/pairplot_iptms2.png")
  """
  return df


def compute_cofolding_metrics(df):
  """
  compute models mean for cofolding empty and plugged ipTM
  """
  df['cofold_mean_empty_i_pTM'] = df[['1_wrepred_empty_i_pTM', '2_wrepred_empty_i_pTM']].mean(axis=1)
  df['cofold_mean_plugged_i_pTM'] = df[['1_wrepred_plugged_i_pTM', '2_wrepred_plugged_i_pTM']].mean(axis=1)
  df['cofold_mean_plugged_Binder_RMSD_to_binding_site'] = df[['1_wrepred_plugged_Binder_RMSD_to_binding_site', '2_wrepred_plugged_Binder_RMSD_to_binding_site']].mean(axis=1)
  df['cofold_mean_empty_Binder_RMSD_to_binding_site'] = df[['1_wrepred_empty_Binder_RMSD_to_binding_site', '2_wrepred_empty_Binder_RMSD_to_binding_site']].mean(axis=1)
  return df


# and a function to extract relevant binders from each final_stat.csv df from bindcraft
def extract_filtered_binders(filtered_df, accepted_folder, filtered_binders_pdbs):
  """ 
  copies the interesting binders in filtered_df from accepted_folder to filtered_binders_pdbs
  """
  print(f"Filters on specific reprediction scores gave {len(filtered_df['Design'])} binders: {filtered_df['Design']}")

  #accepted_folder = f"{design_folder}/Accepted"
  str_binders=''
  for binder_name in filtered_df['Design']:
  # retrieve .pdb file using glob to find the model number
   pdb_files = glob.glob(f"{accepted_folder}/{binder_name}_model*.pdb")
   if not pdb_files:
       print(f"[Warning] No PDB file found for {binder_name} in {accepted_folder}")
       continue # Skip to the next binder if no file is found
   pdb_path = pdb_files[0] # Assuming there's only one matching file, take the first one
  #print(pdb_path)
   str_binders+=f"{pdb_path} "
  # also copy the relevant pdbs to a filtered_binders folder:
   shutil.copy(pdb_path, filtered_binders_pdbs)

  print(str_binders)







