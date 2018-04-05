#!/bin/bash

# Run snakemake
snakemake  -j --snakefile Snakefile.py \
  --cluster-config cluster.json

mv __init__.py scripts/

# Dry run
snakemake -np -j --snakefile Snakefile.py \
  --cluster "sbatch -p testing -t 02:00:00"

# Create graph
snakemake --dag -j --snakefile Snakefile.py \
  --cluster "sbatch -p testing -t 00:30:00" | dot -Tsvg > graph/dag.svg

# Delete all files
rm $(snakemake --snakefile Snakefile.py --summary | tail -n+2 | cut -f1)

sacctmgr show association where user=taavi74

git remote set-url origin git@github.com:avilab/vs.git


conda install -c bioconda snakemake
conda install -c bioconda fastp
conda install -c bioconda fastq-join
conda install -c bioconda cd-hit
conda install -c bioconda tantan

conda install -c bioconda pyfasta
grep '>' output/tantan_goodreads/I1164_12629_Harvard_SIV_196_06_2_24_12_mini.tantan.goodseq.fa | wc -l
conda install -c bioconda ucsc-fasplit


faSplit sequence  {params.n} {output.stub}

# get a password:
wget http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz

conda env export > envs/environment.yml
