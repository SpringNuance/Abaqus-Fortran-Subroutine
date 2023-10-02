#!/bin/bash -l
# Author: Xuan Binh
#SBATCH --job-name=abaqus_subroutine
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --partition=test
#SBATCH --account=project_2007935
#SBATCH --mail-type=ALL
#SBATCH --mail-user=binh.nguyen@aalto.fi

# This script runs in parallel Abaqus example e1 on Puhti server using 10 cores.

unset SLURM_GTIDS

module load abaqus/2022

# Old Intel compilers
module load intel-oneapi-compilers-classic
# module load gcc

cd $PWD

CPUS_TOTAL=$(( $SLURM_NTASKS*$SLURM_CPUS_PER_TASK ))

mkdir tmp_$SLURM_JOB_ID

abq2022 job=geometry input=geometry.inp user=UMAT cpus=$CPUS_TOTAL -verbose 2 standard_parallel=all scratch=tmp_$SLURM_JOB_ID interactive