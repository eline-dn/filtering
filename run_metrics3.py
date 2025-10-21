"""
The purpose of this script is to apply a set of metrics/filters to all the binders in a folder (.pdb format) in order to extract the most relevant ones for experimental validation
1st argument should be the path to the BindCraft output folder, with a final_design_stats.csv and an "Accepted" folder.
2nd argument should be the path to an output folder, were all the analyses can be stored

"""
from rosetta_functions import *
import os
from metrics_utils import *
import pandas as pd
import glob 
import time
from cofolding_utils import *

script_start_time = time.time()
design_folder=sys.argv[1]
output_folder=sys.argv[2]

final_csv = pd.read_csv(f"{design_folder}/final_design_stats.csv")
accepted_folder = f"{design_folder}/Accepted"
params = '/work/lpdi/users/goldbach/software/colabdesign/params'

# empty target:
empty_target_path_orig="/work/lpdi/users/eline/binderdesign/cagedH6_empty.pdb"
# plugged target
plugged_target_path_orig="/work/lpdi/users/eline/binderdesign/trimmed.pdb"


required_cols = [
        "Design", "Sequence","InterfaceResidues","Average_pLDDT", "Average_pTM", "Average_i_pTM",
        "Average_pAE", "Average_i_pAE", "Average_i_pLDDT", "Average_ss_pLDDT",
        "Average_dSASA", "Average_Interface_SASA_%", "Average_Interface_Hydrophobicity",
        "Average_n_InterfaceResidues", "Average_Binder_pLDDT", "Average_Binder_pTM",
        "Average_Binder_pAE", "Average_Binder_RMSD"
    ]

rosetta_cols=['binder_score',
  'surface_hydrophobicity',
  'interface_sc',
    'interface_packstat',
    'interface_dG',
    'interface_dSASA',
    'interface_dG_SASA_ratio',
    'interface_fraction',
    'interface_hydrophobicity',
    'interface_nres',
    'interface_interface_hbonds',
    'interface_hbond_percentage',
              'design_interface_AA_count_types',
              'design_interface_residues']

specific_reprediction_cols=['specific_1_empty_pLDDT', 'specific_1_empty_pTM', 'specific_1_empty_i_pTM', 'specific_1_empty_pAE',
       'specific_1_empty_i_pAE', 'specific_1_empty_Binder_RMSD_to_binding_site', 'specific_2_empty_pLDDT',
       'specific_2_empty_pTM', 'specific_2_empty_i_pTM', 'specific_2_empty_pAE', 'specific_2_empty_i_pAE',
       'specific_2_empty_Binder_RMSD_to_binding_site', 'specific_1_plugged_pLDDT',
       'specific_1_plugged_pTM', 'specific_1_plugged_i_pTM', 'specific_1_plugged_pAE', 'specific_1_plugged_i_pAE',
       'specific_1_plugged_Binder_RMSD_to_binding_site', 'specific_2_plugged_pLDDT',
       'specific_2_plugged_pTM', 'specific_2_plugged_i_pTM', 'specific_2_plugged_pAE', 'specific_2_plugged_i_pAE',
       'specific_2_plugged_Binder_RMSD_to_binding_site']
reprediction_cols=['1_empty_pLDDT', '1_empty_pTM', '1_empty_i_pTM', '1_empty_pAE',
       '1_empty_i_pAE', '1_empty_Binder_RMSD_to_binding_site', '2_empty_pLDDT',
       '2_empty_pTM', '2_empty_i_pTM', '2_empty_pAE', '2_empty_i_pAE',
       '2_empty_Binder_RMSD_to_binding_site', '1_plugged_pLDDT',
       '1_plugged_pTM', '1_plugged_i_pTM', '1_plugged_pAE', '1_plugged_i_pAE',
       '1_plugged_Binder_RMSD_to_binding_site', '2_plugged_pLDDT',
       '2_plugged_pTM', '2_plugged_i_pTM', '2_plugged_pAE', '2_plugged_i_pAE',
       '2_plugged_Binder_RMSD_to_binding_site']


whole_reprediction_cols=['1_wrepred_empty_pLDDT', '1_wrepred_empty_pTM',
       '1_wrepred_empty_i_pTM', '1_wrepred_empty_pAE',
       '1_wrepred_empty_i_pAE','1_wrepred_empty_Binder_RMSD_to_binding_site', '2_wrepred_empty_pLDDT',
       '2_wrepred_empty_pTM', '2_wrepred_empty_i_pTM',
       '2_wrepred_empty_pAE', '2_wrepred_empty_i_pAE','2_wrepred_empty_Binder_RMSD_to_binding_site',
    '1_wrepred_plugged_pLDDT', '1_wrepred_plugged_pTM',
       '1_wrepred_plugged_i_pTM', '1_wrepred_plugged_pAE',
       '1_wrepred_plugged_i_pAE','1_wrepred_plugged_Binder_RMSD_to_binding_site', '2_wrepred_plugged_pLDDT',
       '2_wrepred_plugged_pTM', '2_wrepred_plugged_i_pTM',
       '2_wrepred_plugged_pAE', '2_wrepred_plugged_i_pAE','2_wrepred_plugged_Binder_RMSD_to_binding_site']



df_metrics=pd.DataFrame(columns= required_cols + rosetta_cols + specific_reprediction_cols +reprediction_cols + whole_reprediction_cols)
csv_file= os.path.join(output_folder, 'metrics.csv')
df_metrics.to_csv(csv_file, index=False)

for binder_name in final_csv['Design']: 
#for binder_name in ["n3_l167_s137405_mpnn1","n3_l101_s821023_mpnn2"]: # used for debugging
  # retrieve .pdb file using glob to find the model number
  pdb_files = glob.glob(f"{accepted_folder}/{binder_name}_model*.pdb")
  if not pdb_files:
      print(f"[Warning] No PDB file found for {binder_name} in {accepted_folder}")
      continue # Skip to the next binder if no file is found
  pdb_path = pdb_files[0] # Assuming there's only one matching file, take the first one


### retrieve interesting BC computations and Rosetta metrics
  BC_metrics_dict= load_bindcraft_metrics(final_csv, binder_name, required_cols)
  # compute rosetta metrics
  design_interface_scores_dict = score_interface(pdb_path, binder_chain="B")

  # add them to the metrics df
  BC_metrics_df = pd.DataFrame(BC_metrics_dict,index=[0])
  rosetta_metrics_df = pd.DataFrame(design_interface_scores_dict,index=[0])
  df= pd.concat([BC_metrics_df, rosetta_metrics_df], axis=1)

### 1) compute ipTM for empty barrel and also see if we can repredict the plugged barrel interaction
  binder_sequence=final_csv.loc[final_csv['Design'] == binder_name, 'Sequence'].iloc[0]
  length=len(binder_sequence)

  prediction_model=compile_prediction_models(hardtarget_mode=False,data_dir=params)
  # Generate specific templates 
  empty_target_path_specific=extract_template_path(pdb_path, empty=True, hardtarget_mode=False, design_name=binder_name,output_folder=output_folder)
  plugged_target_path_specific=extract_template_path(pdb_path, empty=False, hardtarget_mode=False, design_name=binder_name, output_folder=output_folder)


  # Run re-predictions with the specific templates
  specific_empty_prediction_stats_df=run_prediction_with_template(model=prediction_model,
                                                         template=empty_target_path_specific, # Use the specific template path
                                                         binder_len=length,
                                                         hardtarget_mode=False,
                                                         binder_sequence=binder_sequence,
                                                         output_folder=output_folder,
                                                         empty=True,
                                                         BC_complex_pdb=pdb_path,
                                                         binder_name=binder_name)
  specific_plugged_prediction_stats_df=run_prediction_with_template(model=prediction_model,
                                                           template=plugged_target_path_specific, # Use the specific template path
                                                           binder_len=length,
                                                           hardtarget_mode=False,
                                                           binder_sequence=binder_sequence,
                                                           output_folder=output_folder,
                                                           empty=False,
                                                           BC_complex_pdb=pdb_path,
                                                           binder_name=binder_name)


  specific_ipTM_reprediction_df=pd.concat([specific_empty_prediction_stats_df, specific_plugged_prediction_stats_df], axis=1)
  df= pd.concat([df, specific_ipTM_reprediction_df], axis=1)
  print(df)
  # if plugged repredicted binder is to far from the original binding site (check RMSD), 
  # run run_prediction_with_template in hardtarget mode
  #ht_prediction_model=compile_prediction_models(hardtarget_mode=True) # predictions will be run using the target + binder complex as a template, 
  #to try if binding results from the BC pipeline cannot be repredicted

### 2) repredict binder structure using the same target template for every structure
  # empty target:
  empty_target_path_orig="/work/lpdi/users/eline/binderdesign/cagedH6_empty.pdb"
  # plugged target
  plugged_target_path_orig="/work/lpdi/users/eline/binderdesign/trimmed.pdb"

  prediction_model=compile_prediction_models(hardtarget_mode=False,data_dir=params) 
  # Run re-predictions with the  non-specific templates
  empty_prediction_stats_df=run_prediction_with_template(model=prediction_model,
                                                         template=empty_target_path_orig, 
                                                         binder_len=length,
                                                         hardtarget_mode=False,
                                                         binder_sequence=binder_sequence,
                                                         output_folder=output_folder,
                                                         empty=True,
                                                         BC_complex_pdb=pdb_path,
                                                         binder_name=binder_name,
                                                         specific=False)
  plugged_prediction_stats_df=run_prediction_with_template(model=prediction_model,
                                                           template=plugged_target_path_orig,
                                                           binder_len=length,
                                                           hardtarget_mode=False,
                                                           binder_sequence=binder_sequence,
                                                           output_folder=output_folder,
                                                           empty=False,
                                                           BC_complex_pdb=pdb_path,
                                                           binder_name=binder_name,
                                                           specific=False)


  ipTM_reprediction_df=pd.concat([empty_prediction_stats_df, plugged_prediction_stats_df], axis=1)
  print(ipTM_reprediction_df)
  df= pd.concat([df, ipTM_reprediction_df], axis=1)
  

### 3) binder and target co-folding, based only on their sequences
  

  empty_reprediction_stats_df=run_whole_reprediction(binder_name=binder_name, binder_sequence=binder_sequence, target_path=plugged_target_path_orig, empty=True, output_folder=output_folder,BC_complex_pdb=pdb_path,params=params)
  plugged_prediction_stats_df=run_whole_reprediction(binder_name=binder_name, binder_sequence=binder_sequence, target_path=empty_target_path_orig, empty=False, output_folder=output_folder,BC_complex_pdb=pdb_path,params=params)
  
  whole_reprediction_df=pd.concat([empty_prediction_stats_df, plugged_prediction_stats_df], axis=1)
  print(whole_reprediction_df)
  df= pd.concat([df, whole_reprediction_df], axis=1)

### eventually: add binder specific data to global metrics df

  df.to_csv(csv_file, mode='a', header=False, index=False)# directly writes the new data in the csv file



# end of the script: how long?
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
print("Finished all designs. Script execution took: "+elapsed_text)
