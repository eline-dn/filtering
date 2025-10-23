#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition=h100
#SBATCH --qos=normal
#SBATCH --gres=gpu:1
#SBATCH --mem 42gb
#SBATCH --time 03:00:00
#SBATCH --output=metrics_%A.log

#SBATCH --mail-type=ALL

# Initialise environment and modules
#source /work/lpdi/users/mpacesa/Pipelines/miniforge3/bin/activate /work/lpdi/users/mpacesa/Pipelines/miniforge3/envs/FastBC
source /work/lpdi/users/mpacesa/Pipelines/miniforge3/bin/activate /work/lpdi/users/mpacesa/Pipelines/miniforge3/envs/BindCraft_kuma
module load gcc/13.2
module load cuda/12.4.1
module load cudnn/8.9.7.29-12

BC_output="$1"
OUTPUT="$2"

python tested_vs_run_metrics.py "${BC_output}" "${OUTPUT}"

