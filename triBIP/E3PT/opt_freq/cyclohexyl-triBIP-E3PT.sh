#!/bin/bash
#SBATCH --job-name=cyclohexyl-triBIP-E3PT
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=pi_hammes_schiffer,day,week
#SBATCH -t 24:00:00

module load Apps/Gaussian/2016-A03

g16 < cyclohexyl-triBIP-E3PT.com > cyclohexyl-triBIP-E3PT.log

rm slurm-${SLURM_JOB_ID}.out
rm core.*
