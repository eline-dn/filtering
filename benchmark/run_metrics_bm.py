"""
goal: apply the metrics to a set of experimentally validated binders to check the filters' predictive power and/ or check reproducibility of the scores
"""
### dependancies
from rosetta_functions import *
import os
from metrics_utils import *
import pandas as pd
import glob 
import time
from cofolding_utils import *
import argparse
from Bio.PDB import PDBList
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer, MMCIFParser


### setup
script_start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("--binder_csv", required =True , type=str, help="path to a csv file with at least the following columns: id (unique for each binder), sequence and target_pdb_id")
parser.add_argument("--output_folder", required=True, type=str, help="to store target templates and analyses")
args = parser.parse_args()

binder_csv = pd.read_csv(args.binder_csv)
output_folder=args.output_folder
os.makedirs(f"{output_folder}/predicted_models", exist_ok=True)
predicted_folder=f"{output_folder}/predicted_models"
params = '/work/lpdi/users/goldbach/software/colabdesign/params'

# ---------------helpers--------------------
# Helper: compute aligned rmsd with superimposer from ca list obtained with get_ca_atoms
def get_ca_atoms(chain):
    atoms = []
    for res in chain:
        if is_aa(res, standard=True) and "CA" in res:
            atoms.append(res["CA"])
    return atoms
    
def sup_rmsd(ref_ca, mov_ca): # compute aligned rmsd with superimposer from ca list obtained with get_ca_atoms
  # ensure same length by trimming the longer one
  L = min(len(ref_ca), len(mov_ca))
  ref_ca = ref_ca[:L]
  mov_ca = mov_ca[:L]

  sup = Superimposer()
  sup.set_atoms(ref_ca, mov_ca)
  rot, tran = sup.rotran
  rmsd = sup.rms
  return(rmsd)
  
## /!\ only for amino acid chains!   
def aligned_chain_rmsd(ref, mov, ref_chain_id, mov_chain_id): # compute one aligned chain rmsd from structures
  first_ref_model = next(ref.get_models())
  first_ref_model_id = first_ref_model.id
  first_mov_model = next(mov.get_models())
  first_mov_model_id = first_mov_model.id
  ref_chain=ref[first_ref_model_id][ref_chain_id]
  mov_chain=mov[first_mov_model_id][mov_chain_id]
  ref_ca = get_ca_atoms(ref_chain)
  mov_ca = get_ca_atoms(mov_chain)
  rmsd=sup_rmsd(ref_ca, mov_ca)
  return rmsd


# fetch target templates
# if not already in output_folder, fetch from pdb
def get_target_template_path(target_name, output_folder):
    path=os.path.join(output_folder, f"{target_name}.pdb")
    if not os.path.isfile(path):
        pdbl = PDBList()
        file_path = pdbl.retrieve_pdb_file(target_name, pdir=output_folder, file_format="pdb")
        print(f"PDB saved to: {file_path}")
        new_path = os.path.join(output_folder, f"{target_name}.pdb")
        shutil.move(file_path, new_path)
        print(f"Renamed to: {new_path}")
        return(new_path)
    else:
        return(path)


# -----------------------iterate on binders --------------------------------
for binder_name in binder_csv['id']: 
    # find target template path, if it doesn't exist fetch it from pbd
    target_pdb_path=get_target_template_path(binder_csv.loc[binder_csv['id'] == binder_name, 'target_pdb_id'].iloc[0],  output_folder)
    
    # 1 - reprediction with pdb template
    binder_sequence=binder_csv.loc[binder_csv['Design'] == binder_name, 'sequence'].iloc[0]
    binder_length=len(binder_sequence)
    prediction_model=compile_prediction_models(hardtarget_mode=False,data_dir=params) 
    # Run re-predictions with the  non-specific templates
    prediction_model.prep_inputs(pdb_filename=target_pdb_path,
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

        predicted_complex_pdb = os.path.join(predicted_folder, f"{binder_name}_model_{model_num+1}_repredicted.pdb")
        prediction_model.save_pdb(predicted_complex_pdb)
        prediction_metrics = copy_dict(prediction_model.aux["log"]) # contains plddt, ptm, i_ptm, pae, i_pae
        # extract the statistics for the model
        stats = {
            f"pLDDT_{model_num}": round(prediction_metrics['plddt'], 2),
            f"pTM_{model_num}": round(prediction_metrics['ptm'], 2),
            f"i_pTM_{model_num}": round(prediction_metrics['i_ptm'], 2),
            f"pAE_{model_num}": round(prediction_metrics['pae'], 2),
            f"i_pAE_{model_num}": round(prediction_metrics['i_pae'], 2)
        }
        # 3 - compute rosetta metrics
        design_interface_scores_dict = score_interface(predicted_complex_pdb, binder_chain="B")
        for key in design_interface_scores_dict.keys():
            key = key + f"_{model_num}"
        # 4 -  compute RMSDs: target rmsd and  later coherence between model 1 and model 2
        rmsds={}
        from Bio.PDB import PDBParser, Selection
        parser = PDBParser(QUIET=True)
        ref=parser.get_structure("a",target_pdb_path)
        mov = parser.get_structure("x", predicted_complex_pdb)
        rmsds[f"target_rmsd_{model_num}"]=aligned_chain_rmsd(ref, mov, ref_chain_id="A", mov_chain_id="A")
        rmsds[f"binder_rmsd_{model_num}"]=aligned_chain_rmsd(ref, mov, ref_chain_id="B", mov_chain_id="B")
        data={**rmsds, **design_interface_scores_dict, **stats}
        prediction_stats = {**prediction_stats, **data}
 
    # convert to df
    prediction_stats["id"]=binder_name
    df=pd.DataFrame(data=prediction_stats, index=[prediction_stats["id"]])
    # mean or compare the two model's outputs:
    
    # save to csv:
    csv_path=f"{outdir}/{prefix}_Colab_Design_reprediction_metrics.csv"
    df.to_csv(csv_path, mode="a", index=False, header=not pd.io.common.file_exists(csv_path))


# end of the script: how long?
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
print("Finished all designs. Script execution took: "+elapsed_text)
