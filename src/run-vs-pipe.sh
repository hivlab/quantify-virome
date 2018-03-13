#!/bin/bash

#The job should run on the testing partition
#SBATCH -p testing

#The name of the job
#SBATCH -J parallel_snakemake

#Number of compute nodes the job requires
#SBATCH -N 2

#The job requires 1 task per node
#SBATCH --ntasks-per-node=8

#The maximum walltime of the job is a half hour
#SBATCH -t 00:30:00

# Send remainders
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tapa741@gmail.com

#SBATCH -o /gpfs/hpchome/taavi74/Projects/vs/output/logs/stdout.%j

#Here we call srun to launch the uname command in parallel
module purge
module load adapterremoval/2.1.7 python-3.6.0 fastq-join

cd /gpfs/hpchome/taavi74/Projects/vs

SAMPLE=I1164_12629_Harvard_SIV_196_06_2_24_12

source env/bin/activate
snakemake --snakefile src/snakefile.py output/stitched_reads/$SAMPLE.stitch-length-report
deactivate

