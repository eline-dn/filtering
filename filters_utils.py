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


