# Filtering: *Selecting binders generated with BindCraft for experimental testing*
## To compute the different reprediction and rosetta metrics:
run 'sbatch run_metrics.slurm' or 'python run_metrics3.py' for each condition
with the following arguments: 
- 1st argument should be the path to the BindCraft output folder, with a final_design_stats.csv and an "Accepted" folder.
- 2nd argument should be the path to an output folder, were all the analyses can be stored

## To filter the binders based on these metrics:
run 'python run_filters.py' with the following arguments:
- path to the folder where every condition's BC output folder can be found. Each condition's folder must contain an Accepted folder with the binders' pdbs
- path to the folder where every condition's run_metrics.py output folder can be found. Each folder has a ./metrics.csv
- where to store the successful binders and plots

/!\ run_metrics needs to be run individually for each condition, while run_filters is a batch run for all the conditions to extract binders from.

