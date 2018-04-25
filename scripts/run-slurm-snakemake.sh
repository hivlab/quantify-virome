#!/bin/bash

# Run snakemake
snakemake -j --cluster-config cluster.json \
             --cluster "sbatch -J {cluster.name} \
             -p {cluster.partition} \
             -t {cluster.time} \
             --mem {cluster.mem} \
             --output {cluster.output}"

# Dry run
snakemake -np -j --cluster-config cluster.json \
             --cluster "sbatch -J {cluster.name} \
             -p {cluster.partition} \
             -t {cluster.time} \
             --mem {cluster.mem} \
             --output {cluster.output}"

# Create graph
snakemake --dag -j --cluster-config cluster.json \
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

git add --all
git commit -m "not necessary"
git push
git config --global push.default simple
git pull

rm scripts/__init__.py
rm output/repeatmasker_good/*.*

source activate virome
pwd > ~/Projects/vs/db.txt
source deactivate
cd ~/Projects/vs
touch samples.csv
mkdir rules
touch rules/munge.smk
touch rules/mask.smk

conda install -c bioconda bwa
conda install -c bioconda picard
export GIRUSER=taavipall >> ~/.bashrc
export GIRPASS=lk563m >> ~/.bashrc
wget --user taavipall --password lk563m http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz

cd /gpfs/software/VirusSeeker/databases/ref_genomes/
cat .snakemake/log/2018-04-10T234450.700272.snakemake.log
sbatch test.sh
squeue -u taavi74
which bwa
