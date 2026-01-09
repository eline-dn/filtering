#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --partition=l40s
#SBATCH --qos=normal
#SBATCH --gres=gpu:1
#SBATCH --mem 30gb
#SBATCH --time 06:00:00
#SBATCH --output=metrics_BC%A.log

#SBATCH --mail-type=ALL

# Initialise environment and modules
#source /work/lpdi/users/mpacesa/Pipelines/miniforge3/bin/activate /work/lpdi/users/mpacesa/Pipelines/miniforge3/envs/FastBC
source /work/lpdi/users/mpacesa/Pipelines/miniforge3/bin/activate /work/lpdi/users/mpacesa/Pipelines/miniforge3/envs/BindCraft_kuma
module load gcc/13.2
module load cuda/12.4.1
module load cudnn/8.9.7.29-12

python benchmark/run_metrics_bm.py --binder_csv "/work/lpdi/users/eline/projects/filtering/benchmark/BC_reprod.csv" --output_folder /work/lpdi/users/eline/projects/filtering/benchmark
