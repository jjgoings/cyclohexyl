#!/bin/bash
#SBATCH --job-name=cyclohexyl-monoBIP-E0PT
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=pi_hammes_schiffer,day,week
#SBATCH -t 08:00:00

module load Apps/Gaussian/2016-A03

g16 < cyclohexyl-monoBIP-E0PT.com > cyclohexyl-monoBIP-E0PT.log

rm slurm-${SLURM_JOB_ID}.out
rm core.*