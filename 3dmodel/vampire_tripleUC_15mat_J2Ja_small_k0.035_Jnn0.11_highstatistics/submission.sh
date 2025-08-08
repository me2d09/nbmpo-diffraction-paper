#!/bin/bash
#SBATCH --no-requeue
#SBATCH --job-name="J2Ja_small_k0.035_Jnn0.11_Ba_HS"
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=490000

source /opt/intel/oneapi/setvars.sh
source /mnt/home/dossan_f/venv/main/bin/activate
PATH=$PATH:$HOME/codes/vampire

python launch.py
