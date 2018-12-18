#!/bin/bash
#SBATCH --job-name=cyclohexyl-diBIP-E2PT
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=pi_hammes_schiffer,day,week
#SBATCH -t 08:00:00

module load Apps/Gaussian/2016-A03

g16 < cyclohexyl-diBIP-E2PT.com > cyclohexyl-diBIP-E2PT.log

rm slurm-${SLURM_JOB_ID}.out
rm core.*