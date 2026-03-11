#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --partition=h100
#SBATCH --qos=normal
#SBATCH --gres=gpu:1
#SBATCH --mem 30gb
#SBATCH --time 72:00:00
#SBATCH --output=metrics_pb%A.log

#SBATCH --mail-type=ALL

# Initialise environment and modules
#source /work/lpdi/users/mpacesa/Pipelines/miniforge3/bin/activate /work/lpdi/users/mpacesa/Pipelines/miniforge3/envs/FastBC
source /work/lpdi/users/mpacesa/Pipelines/miniforge3/bin/activate /work/lpdi/users/mpacesa/Pipelines/miniforge3/envs/BindCraft_kuma
module load gcc/13.2
module load cuda/12.4.1
module load cudnn/8.9.7.29-12

python /work/lpdi/users/eline/projects/filtering/benchmark/run_bm_pb.py --binder_csv "/work/lpdi/users/eline/projects/filtering/benchmark/protbase/proteinbase_all_data_28_10_2025_cleaned.csv" --output_folder /work/lpdi/users/eline/projects/filtering/benchmark/protbase
