"""
goal: apply the metrics to a set of experimentally validated binders to check the filters' predictive power and/ or check reproducibility of the scores
here: run on the proteinbase database by adaptyv bio
"""
### dependencies
import os, sys
SCRIPT_PATH = os.path.dirname(__file__)
sys.path.append(f"{SCRIPT_PATH}/..")
import shutil
from rosetta_functions import *
from metrics_utils import *
import pandas as pd
import glob 
import time
from cofolding_utils import *
import argparse
from Bio.PDB import PDBList
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer, MMCIFParser
from collections import defaultdict
import numbers

### setup
script_start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("--binder_csv", required =True , type=str, help="path to a csv file with at least the following columns: id (unique for each binder), sequence and binding_target")
parser.add_argument("--output_folder", required=True, type=str, help="to store target templates and analyses")
args = parser.parse_args()

binder_csv = pd.read_csv(args.binder_csv)
output_folder=args.output_folder
os.makedirs(f"{output_folder}/predicted_models", exist_ok=True)
predicted_folder=f"{output_folder}/predicted_models"
params = '/work/lpdi/users/goldbach/software/colabdesign/params'

# prepare the input csv by adding the target_pdb_id column:
# initial target pdb file associated with each binder:
target_dict={# from bc:
  "IFNAR2": "2LAG",
            "CD45" : "5FMV", 
            "BBF-14": "9HAG",
            "DerF7": "3UV1",
            "DerF21":"5YNY",
              "CbAgo": "6QZK",
            "HER2":"1N8Z",
            "SpCas9": "4ZT0",
            # from proteinbase
  "human-pdgfr-beta":"3MJG", #ok
  "human-insulin-receptor":"9DNI", #ok
  "human-gm2a":"1PUB", #ok
  "human-orm2": "3KQ0", #ok
  "human-ambp": "3QKG",
  "spcas9": "4ZT0",
  "der21": "5YNY",
  "der7": "3UV1",
  "ifnar2": "2LAG",
"human-mzb1-perp1":"7AAH", #ok
"il7r": "5J11", #ok
"human-idi2":"2PNY", #ok
"hnmt":"1JQD", #ok
"human-pmvk":"3CH4", # failed
"human-phyh":"2A1X",
"egfr":"1IVO"
}
target2skip= ["BetV1" , "PD1", "PD-L1","CLDN1","Sas6", "human-tnfa", "mdm2", "fgf-r1"]

def target2pdb(target_name, target_dict, target2skip):
  if (target_name in target2skip) or (target_name not in target_dict.keys()) :
    return(np.nan)
  else:
    return target_dict[target_name]
  
binder_csv["target_pdb_id"]=binder_csv["binding_target"].apply(target2pdb, target_dict=target_dict, target2skip=target2skip)
binder_csv=binder_csv.dropna(subset=["target_pdb_id"])
binder_csv.to_csv(f"{output_folder}/ref_csv.csv")

# ---------------helpers--------------------
# relaxation
# clean unnecessary rosetta information from PDB
def clean_pdb(pdb_file):
    # Read the pdb file and filter relevant lines
    with open(pdb_file, 'r') as f_in:
        relevant_lines = [line for line in f_in if line.startswith(('ATOM', 'HETATM', 'MODEL', 'TER', 'END', 'LINK'))]

    # Write the cleaned lines back to the original pdb file
    with open(pdb_file, 'w') as f_out:
        f_out.writelines(relevant_lines)
# Relax designed structure
def pr_relax(pdb_file, relaxed_pdb_path):
    if not os.path.exists(relaxed_pdb_path):
        # Generate pose
        pose = pr.pose_from_pdb(pdb_file)
        start_pose = pose.clone()

        ### Generate movemaps
        mmf = MoveMap()
        mmf.set_chi(True) # enable sidechain movement
        mmf.set_bb(True) # enable backbone movement, can be disabled to increase speed by 30% but makes metrics look worse on average
        mmf.set_jump(False) # disable whole chain movement

        # Run FastRelax
        fastrelax = FastRelax()
        scorefxn = pr.get_fa_scorefxn()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf) # set MoveMap
        fastrelax.max_iter(200) # default iterations is 2500
        fastrelax.min_type("lbfgs_armijo_nonmonotone")
        fastrelax.constrain_relax_to_start_coords(True)
        fastrelax.apply(pose)

        # Align relaxed structure to original trajectory
        align = AlignChainMover()
        align.source_chain(0)
        align.target_chain(0)
        align.pose(start_pose)
        align.apply(pose)

        # Copy B factors from start_pose to pose
        for resid in range(1, pose.total_residue() + 1):
            if pose.residue(resid).is_protein():
                # Get the B factor of the first heavy atom in the residue
                bfactor = start_pose.pdb_info().bfactor(resid, 1)
                for atom_id in range(1, pose.residue(resid).natoms() + 1):
                    pose.pdb_info().bfactor(resid, atom_id, bfactor)

        # output relaxed and aligned PDB
        pose.dump_pdb(relaxed_pdb_path)
        clean_pdb(relaxed_pdb_path)
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


def align_to_chain(ref,mov,mapping): # align a structure to a reference structure based on the chains in mapping
  # concat all the chains' CA in the same order
  # provide in this case a mapping of the ref and mov structure chains in the shape
  # { ref_chain_id1: mov_chain_id1, ref_chain_id2: mov_chain_id2,...}
  # !!! the targeted chains should be the same lentgh!!
  ref_ca_list=list()
  mov_ca_list=list()
  first_ref_model = next(ref.get_models())
  first_ref_model_id = first_ref_model.id
  first_mov_model = next(mov.get_models())
  first_mov_model_id = first_mov_model.id
  for ref_chain_id, mov_chain_id in mapping.items():
    ref_chain=ref[first_ref_model_id][ref_chain_id]
    mov_chain=mov[first_mov_model_id][mov_chain_id]
    ref_ca = get_ca_atoms(ref_chain)
    mov_ca = get_ca_atoms(mov_chain)
    ref_ca_list+=ref_ca
    mov_ca_list+=mov_ca

  if len(ref_ca_list) != len(mov_ca_list):
    raise ValueError(
        f" the list of atoms to align are of different lengths, check the chain mapping and your chain lengths"
        f"(number of ref atoms {len(ref_ca_list)} and number of mov atoms {len(mov_ca_list)} )."
    )
  sup = Superimposer()
  sup.set_atoms(ref_ca_list, mov_ca_list)
  rot, tran = sup.rotran
  # Apply transform to ALL atoms in the moving structure
  for atom in mov.get_atoms():
      atom.transform(rot, tran)
  # return the modified mov structure: (the structure will be modified anyway, even if you use a copy after that
  return(mov)

def unaligned_rmsd(ref, mov, mapping): # compute unaligned rmsd for the chains in mapping, not for ligands!
  def ca_map(chain):
    out = {}
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        if "CA" in res:
            hetflag, resseq, icode = res.get_id()
            out[(resseq, icode)] = res["CA"]
    return out

  first_ref_model = next(ref.get_models())
  first_ref_model_id = first_ref_model.id
  first_mov_model = next(mov.get_models())
  first_mov_model_id = first_mov_model.id
  ref_ca_list={}
  mov_ca_list={}
  for ref_chain_id, mov_chain_id in mapping.items():
    ref_chain=ref[first_ref_model_id][ref_chain_id]
    mov_chain=mov[first_mov_model_id][mov_chain_id]
    ref_ca = ca_map(ref_chain)
    mov_ca = ca_map(mov_chain)
    ref_ca_list={**ref_ca_list, **ref_ca}
    mov_ca_list={**mov_ca_list, **mov_ca}
  
  # Common residue keys, ordered by residue number then insertion code
  common = sorted(set(ref_ca_list.keys()).intersection(mov_ca_list.keys()),
                  key=lambda k: (k[0], (k[1] or " ")))
  
  if len(common) < 3:
      raise ValueError(
          f"Not enough matched residues with Cα to compute RMSD without alignment "
          f"(found {len(common)})."
      )
  
  ref_coords = np.array([ref_ca_list[k].get_coord() for k in common], dtype=float)
  mov_coords = np.array([mov_ca_list[k].get_coord() for k in common], dtype=float)
  
  # Unaligned RMSD = sqrt(mean(||ref - mov||^2))
  diffs = ref_coords - mov_coords
  rmsd = float(np.sqrt((diffs * diffs).sum(axis=1).mean()))
  print(round(rmsd, 2))
  return(rmsd)


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


from collections import defaultdict

def average_paired_metrics(d):
    grouped = defaultdict(list)

    # group values by base metric name
    for key, value in d.items():
        if "_" in key and key.rsplit("_", 1)[-1].isdigit():
            base, suffix = key.rsplit("_", 1)
            # keep only numeric values
            if isinstance(value, numbers.Number):
                grouped[base].append(value)

    # compute averages
    averaged = {}
    for base, values in grouped.items():
        if len(values) != 2:
            raise ValueError(
                f"Metric '{base}' does not have exactly two values: {values}"
            )
        averaged[base] = sum(values) / 2

    return averaged

# -----------------------iterate on binders --------------------------------
for binder_name in binder_csv['id']: 
    # find target template path, if it doesn't exist fetch it from pbd
    target_pdb_path=get_target_template_path(binder_csv.loc[binder_csv['id'] == binder_name, 'target_pdb_id'].iloc[0],  output_folder)
  
    # check that chain A exists in  ref pdb:
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    ref=parser.get_structure("a",target_pdb_path)
    ref_model = next(ref.get_models())
    try:
      ref_chain = ref_model["A"]
    except KeyError:
      print(f"Reference chain 'A' not found in {target_pdb_path}.")
      continue
    # check if this binder was already done:
    predicted_complex_pdb = os.path.join(predicted_folder, f"{binder_name}_model_1_repredicted.pdb")
    if os.path.isfile(predicted_complex_pdb):
        print(f"Binder {binder_name} already processed, skipping.")
        continue
    # 1 - reprediction with pdb template
    binder_sequence=binder_csv.loc[binder_csv['id'] == binder_name, 'sequence'].iloc[0]
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
        # add relaxation before computing the metrics
        # Relax binder to calculate statistics
        predicted_complex_relaxed = os.path.join(predicted_folder, f"{binder_name}_model_{model_num+1}_relaxed.pdb")
        pr_relax(predicted_complex_pdb, predicted_complex_relaxed)

        
        # 3 - compute rosetta metrics
        design_interface_scores_dict = score_interface(predicted_complex_relaxed, binder_chain="B")
        design_interface_scores_dict = {
        f"{key}_{model_num}": value
        for key, value in design_interface_scores_dict.items()
        }

        # 4 -  compute RMSDs: target rmsd and  later coherence between model 1 and model 2
        #rmsds={}
        from Bio.PDB import PDBParser, Selection
        parser = PDBParser(QUIET=True)
        ref=parser.get_structure("a",target_pdb_path)
        mov = parser.get_structure("x", predicted_complex_pdb)
        #rmsds[f"target_rmsd_{model_num}"]=aligned_chain_rmsd(ref, mov, ref_chain_id="A", mov_chain_id="A")
        #rmsds[f"binder_rmsd_{model_num}"]=aligned_chain_rmsd(ref, mov, ref_chain_id="B", mov_chain_id="B")
        data={**design_interface_scores_dict, **stats}
        prediction_stats = {**prediction_stats, **data}
        prediction_stats[f"target_rmsd_{model_num}"]=aligned_chain_rmsd(ref, mov, ref_chain_id="A", mov_chain_id="A")
        #print(model_num, prediction_stats)

    # mean or compare the two model's outputs:
    print(prediction_stats)
    avg_metrics = average_paired_metrics(prediction_stats)
    prediction_stats={**prediction_stats, **avg_metrics}

    # compute distance between the two predicted binding sites
    ref=parser.get_structure("a",os.path.join(predicted_folder, f"{binder_name}_model_1_repredicted.pdb"))
    mov = parser.get_structure("x", os.path.join(predicted_folder, f"{binder_name}_model_2_repredicted.pdb"))
    prediction_stats["binder_delta_rmsd"]=aligned_chain_rmsd(ref, mov, ref_chain_id="B", mov_chain_id="B")
    mov_aligned=align_to_chain(ref,mov,{"A":"A", "B":"B"})
    prediction_stats["d_binding_site"]=unaligned_rmsd(ref, mov_aligned, {"B":"B"})

    # convert to df
    prediction_stats["id"]=binder_name
    df=pd.DataFrame(data=prediction_stats, index=[prediction_stats["id"]]) 
    # save to csv:
    csv_path=os.path.join(output_folder, "binder_scoring_pbase.csv")
    df.to_csv(csv_path, mode="a", index=False, header=not pd.io.common.file_exists(csv_path))


# end of the script: how long?
elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
print("Finished all designs. Script execution took: "+elapsed_text)
