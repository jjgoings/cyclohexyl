#!/bin/bash
#SBATCH --job-name=cyclohexyl-monoBIP-TS2_rev
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=pi_hammes_schiffer,day,week
#SBATCH -t 24:00:00

module load Apps/Gaussian/2016-A03

g16 < cyclohexyl-monoBIP-TS2_rev.com > cyclohexyl-monoBIP-TS2_rev.log

rm slurm-${SLURM_JOB_ID}.out
rm core.*
