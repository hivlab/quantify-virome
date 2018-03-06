#!/bin/bash

#The job should run on the testing partition
#SBATCH -p testing

#The name of the job
#SBATCH -J parallel_uname

#The job requires 4 compute nodes
#SBATCH -N 4

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 00:30:00

# Send remainders
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tapa741@gmail.com

#SBATCH -o /gpfs/hpchome/taavi74/Projects/vs/output/stdout.%j
#SBATCH -e /gpfs/hpchome/taavi74/Projects/vs/output/stderr.%j

#Here we call srun to launch the uname command in parallel
srun uname -a
sleep 30

