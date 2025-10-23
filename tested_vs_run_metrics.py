"""
The purpose of this script is to apply a set of metrics/filters to all the binders that have been experimentally tested in the BindCraft paper
1st argument should be the path to a folder, with a design_stats.csv file with binders and binding information and an "Accepted" folder, with the binders' pdbs. 
2nd argument should be the path to an output folder, were all the analyses can be stored

"""
#from rosetta_functions import *
import os, shutil
from metrics_utils import *
import pandas as pd
import glob 
import time
from cofolding_utils import *
from Bio.PDB import PDBList


script_start_time = time.time()
design_folder=sys.argv[1]
output_folder=sys.argv[2]

final_csv = pd.read_csv(f"{design_folder}/BC_design_stats.csv")
accepted_folder = f"{design_folder}/Accepted"
params = '/work/lpdi/users/goldbach/software/colabdesign/params'
# initial target pdb file associated with each binder:
target_dict={"IFNAR2": "2LAG",
            "CD45" : "5FMV", 
            "BBF-14": "9HAG",
            "DerF7": "3UV1",
            "DerF21":"5YNY",
            "CbAgo": "6QZK",
            "HER2":"1N8Z",
            "SpCas9": "4ZT0"}
target2skip= ["BetV1" , "PD1", "PD-L1"]
# BetV1: 7MXL
# create: target_path_dict, (fetch all target pdb structures and save them as a .pdb file
pdbl = PDBList()
target_path_dict={}
for target, pdb_id in target_dict.items():
            target_dir=output_folder
            # download structure (.pdb file)
            file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=target_dir, file_format="pdb")
            print(f"PDB saved to: {file_path}")
            new_path = os.path.join(target_dir, f"{target}.pdb")
            shutil.move(file_path, new_path)
            print(f"Renamed to: {new_path}")
            target_path_dict[target]=new_path

print(target_path_dict)


required_cols = [
        "Target","DesignName","Binding","Sequence","Average_pLDDT", "Average_pTM", "Average_i_pTM",
        "Average_pAE", "Average_i_pAE", "Average_i_pLDDT", "Average_ss_pLDDT",
        "Average_dSASA", "Average_Interface_SASA_%", "Average_Interface_Hydrophobicity",
        "Average_n_InterfaceResidues", "Average_Binder_pLDDT", "Average_Binder_pTM",
        "Average_Binder_pAE", "Average_Binder_RMSD"
    ]



reprediction_cols=['1_it_pLDDT', '1_it_pTM', '1_it_i_pTM', '1_it_pAE',
       '1_it_i_pAE', '1_it_Binder_RMSD_to_binding_site', '2_it_pLDDT',
       '2_it_pTM', '2_it_i_pTM', '2_it_pAE', '2_it_i_pAE',
       '2_it_Binder_RMSD_to_binding_site']


df_metrics=pd.DataFrame(columns= required_cols  + reprediction_cols )
csv_file= os.path.join(output_folder, 'metrics.csv')
df_metrics.to_csv(csv_file, index=False)

for binder_name in final_csv['DesignName']: 
#for binder_name in ["n3_l167_s137405_mpnn1","n3_l101_s821023_mpnn2"]: # used for debugging
  # retrieve .pdb file using glob to find the model number
  pdb_files = glob.glob(f"{accepted_folder}/{binder_name}.pdb")
  if not pdb_files:
      print(f"[Warning] No PDB file found for {binder_name} in {accepted_folder}")
      continue # Skip to the next binder if no file is found
  pdb_path = pdb_files[0] # Assuming there's only one matching file, take the first one

 
### retrieve interesting BC computations and Rosetta metrics
  BC_metrics_dict= load_bindcraft_metrics_bis(final_csv, binder_name, required_cols)
  """    
  # compute rosetta metrics
  design_interface_scores_dict = score_interface(pdb_path, binder_chain="B")

  # add them to the metrics df
  BC_metrics_df = pd.DataFrame(BC_metrics_dict,index=[0])
  rosetta_metrics_df = pd.DataFrame(design_interface_scores_dict,index=[0])
  df= pd.concat([BC_metrics_df, rosetta_metrics_df], axis=1)
  """
  df = pd.DataFrame(BC_metrics_dict,index=[0])

  binder_sequence=final_csv.loc[final_csv['DesignName'] == binder_name, 'Sequence'].iloc[0]
  binder_length=len(binder_sequence)
  # get target pdb name from the target_dictionnary:
  
  target_name=final_csv.loc[final_csv['DesignName'] == binder_name, 'Target'].iloc[0]
  if target_name in target2skip:
          continue # skip to the next binder
  else:
          # fetch pdb path from the target_path_dict:
          initial_target_pdb=target_path_dict[target_name]
                   

### 2) repredict binder structure using the same target template for every structure, i.e. the initial target
  

  prediction_model=compile_prediction_models(hardtarget_mode=False,data_dir=params) 
  # Run re-predictions with the  non-specific templates
  prediction_model.prep_inputs(pdb_filename=initial_target_pdb,
                        chain="A",
                        #binder_chain="B",# do not specifiy if the template only contains the target
                        binder_len=binder_length,
                        rm_target_seq=False, #b
                        use_binder_template=False, #a
                        rm_template_ic=False #c
                        )
  prediction_stats = {}
  for model_num in [0,1]:
     prediction_model.predict(seq=binder_sequence,
                    models=[model_num],
                    num_recycles=3)
     os.makedirs(f"{output_folder}/predicted_models", exist_ok=True)
     predicted_folder=f"{output_folder}/predicted_models"
     predicted_complex_pdb = os.path.join(predicted_folder, f"{binder_name}_model_{model_num+1}_repredicted.pdb")
     prediction_model.save_pdb(predicted_complex_pdb)
     prediction_metrics = copy_dict(prediction_model.aux["log"]) # contains plddt, ptm, i_ptm, pae, i_pae

            # extract the statistics for the model
     stats = {
                f"it_pLDDT": round(prediction_metrics['plddt'], 2),
                f"it_pTM": round(prediction_metrics['ptm'], 2),
                f"it_i_pTM": round(prediction_metrics['i_ptm'], 2),
                f"it_pAE": round(prediction_metrics['pae'], 2),
                f"it_i_pAE": round(prediction_metrics['i_pae'], 2)
            }
     # unaligned RMSD calculate to determine if binder is in the designed binding site
     rmsd_site = unaligned_rmsd(pdb_path, predicted_complex_pdb, "B", "B")
     stats[f"it_Binder_RMSD_to_binding_site"] = rmsd_site # this should be used to filter the models that are binding in the predicted binding site

     prediction_stats[model_num+1] = stats # 2 dictionnaries index 1 and 2 to eventually add to the metrics df
  ipTM_reprediction_df = transform_prediction_stats_to_df(prediction_stats)


  print(ipTM_reprediction_df)
  df= pd.concat([df, ipTM_reprediction_df], axis=1)

  """
  ### 3) binder and target co-folding, based only on their sequences


  empty_reprediction_stats_df=run_whole_reprediction(binder_name=binder_name, binder_sequence=binder_sequence, target_path=empty_target_path_orig, empty=True, output_folder=output_folder,BC_complex_pdb=pdb_path,params=params)
  plugged_prediction_stats_df=run_whole_reprediction(binder_name=binder_name, binder_sequence=binder_sequence, target_path=plugged_target_path_orig, empty=False, output_folder=output_folder,BC_complex_pdb=pdb_path,params=params)

  whole_reprediction_df=pd.concat([empty_prediction_stats_df, plugged_prediction_stats_df], axis=1)
  print(whole_reprediction_df)
  df= pd.concat([df, whole_reprediction_df], axis=1)
  """

  ### eventually: add binder specific data to global metrics df
  df.to_csv(csv_file, mode='a', header=False, index=False)# directly writes the new data in the csv file



# end of the script: how long?
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
print("Finished all designs. Script execution took: "+elapsed_text)
