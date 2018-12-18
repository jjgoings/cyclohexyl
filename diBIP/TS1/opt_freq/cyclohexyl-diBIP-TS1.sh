#!/bin/bash
#SBATCH --job-name=cyclohexyl-diBIP-TS1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=pi_hammes_schiffer,day,week
#SBATCH -t 24:00:00

module load Apps/Gaussian/2016-A03

g16 < cyclohexyl-diBIP-TS1.com > cyclohexyl-diBIP-TS1.log

rm slurm-${SLURM_JOB_ID}.out
rm core.*