#!/bin/bash

source activate virome

# Run snakemake
snakemake -j --use-conda --cluster-config cluster.json  \
             --cluster "sbatch -J {cluster.name} \
             -p {cluster.partition} \
             -t {cluster.time} \
             --mem {cluster.mem} \
             --output {cluster.output} \
             --cpus-per-task {cluster.cpus-per-task}"

# Dry run
snakemake -n

# Create graph
snakemake --dag | dot -Tsvg > graph/dag.svg

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

/gpfs/hpchome/taavi74/Downloads/repbase.info

export REPEATMASKER_REPBASE_FILE="/gpfs/hpchome/taavi74/databases/Repbase/Libraries/RMRBSeqs.embl"
TRF_DIR=$CONDA_PREFIX/bin
RMBLAST_DIR=$CONDA_PREFIX/bin
HMMER_DIR=$CONDA_PREFIX/bin
echo $REPEATMASKER_REPBASE_FILE
echo $RMBLAST_DIR
export RMBLAST_DIR=/gpfs/hpchome/taavi74/miniconda3/envs/snakemake/bin

$CONDA_PREFIX/bin/rmblastn

snakemake -j --use-conda --snakefile Snakemake_repmask

 mv 08_split_fasta/*.fa.* 09_repeatmasker

/gpfs/hpchome/taavi74/miniconda3/envs/vir/bin/


rm ~/miniconda3/envs/virome/share/RepeatMasker/RepeatMaskerConfig.pm


    mkdir -p $CONDA_PREFIX/share/RepeatMasker/
    cp envs/RepeatMaskerConfig.pm $CONDA_PREFIX/share/RepeatMasker/

/gpfs/hpchome/taavi74/Projects/vs/.snakemake/conda/02ce30bd
source deactivate
source activate 02ce30bd
the assessment grid

echo $REPEATMASKER_REPBASE_FILE
-------------------


mv {params.dir}/08_split_fasta/*.fa.* {params.dir}/09_repeatmasker
    cd {params.dir}/09_repeatmasker
    if [ ! -n "$(find . -maxdepth 1 -name 'tantan.goodseq.*.fa.masked' -print -quit)" ]
    then
       touch tantan.goodseq.1.fa.masked
    fi

rm $CONDA_PREFIX/share/RepeatMasker/RepeatMaskerConfig.pm

samples = pd.read_table("samples.tsv", sep = "\s+", index_col = "sample", dtype = str)
print(samples.columns.tolist())
class()

samples.index.get_values()

samples["sample"].tolist()
s = samples["sample"].
samples.iloc[0,0]

cd /gpfs/hpchome/taavi74/databases/Repbase/Libraries/

RepeatMasker -qq -pa 8 output/I1164_12629_Harvard_SIV_196_06_2_24_12_mini/08_split_fasta/tantan.goodseq.1.fa -dir output/I1164_12629_Harvard_SIV_196_06_2_24_12_mini/09_repeatmasker

echo $(zcat data/I1164_12629_Harvard_SIV_196_06_2_24_12_SE2.fastq.gz|wc -l)/4|bc
1236403


conda create -n virome -c bioconda snakemake repeatmasker trf rmblast
conda info --envs
conda env remove -n rstudio


############################################################
# https://git.wageningenur.nl/warri004/makerAnnotation/tree/d761cb62a785352ee8ed052877834493c1f9c912
# conf files excerpt
"executables":
	{
		"rmblast": "{base_dir}{executables}rmblast/bin/rmblastn",
		"makeblastdb": "{base_dir}{executables}blast_plus/bin/makeblastdb",
		"blastn": "{base_dir}{executables}blast_plus/bin/blastn",
		"blastp": "{base_dir}{executables}blast_plus/bin/blastp",
		"blastx": "{base_dir}{executables}blast_plus/bin/blastx",
		"tblastx": "{base_dir}{executables}blast_plus/bin/tblastx",
		"RepeatMasker": "{base_dir}{executables}repeatmasker/RepeatMasker",
		"repeatscout": "{base_dir}{executables}repeatscout/",
		"repeatmodeler": "{base_dir}{executables}repeatmodeler/",
		"exonerate": "{base_dir}{executables}exonerate/bin/exonerate",
		"snap": "{base_dir}{executables}snap/snap",
		"augustus": "{base_dir}{executables}augustus/bin/augustus",
		"trf": "{base_dir}{executables}trf/trf",
		"repbase": "{base_dir}{executables}repbase/",
		"maker": "{base_dir}{executables}maker/bin/maker",
		"perl" : "{home_dir}perl5/perlbrew/perls/perl-5.18.0/bin/perl",
		"recon" : "{base_dir}{executables}recon/bin",
		"cpanm" : "{home_dir}perl5/perlbrew/bin/cpanm"
	}


--------------------------
# Download and install GIRI RepBase
rule repbase_installation:
    output: CONFIG["executables"]["repbase"]
    shell:
        "TEMPDIR=`mktemp -p " + CONFIG["base"]["download_dir"] + " -d` && "
        "cd $TEMPDIR && "
        "curl --remote-name --anyauth --netrc-file "+CONFIG['executable_sources']['repeatmaskerlibraries_URL_netrc']+" "+CONFIG['executable_sources']['repeatmaskerlibraries_URL']+" && "
        "tar xzf repeatmaskerlibraries* && "
        "cp -R Libraries {output} && "
        "cp -R repeatmaskerlibraries* {output} && "
        "cd / && rm -rf $TEMPDIR && "
        "touch {output};"

# Download and install Repeatmasker
rule repeatmasker_installation:
    input:
        #REPBASE=rules.repbase_installation.output,
        RMBLAST=rules.rmblast_installation.output,
        TRF=rules.trf_installation.output
    output: CONFIG["executables"]["RepeatMasker"]
    params: dir = strip_path_level(CONFIG["executables"]["RepeatMasker"],1), RMBLAST = strip_path_level(CONFIG["executables"]["rmblast"],1)
    shell:
        "TEMPDIR=`mktemp -p " + CONFIG["base"]["download_dir"] + " -d` && "
        "cd $TEMPDIR && "
        "wget"+CONFIG['executable_sources']['wget_options']+CONFIG['executable_sources']['repeatmasker_URL']+" && "
        "tar xzf RepeatMasker-open-* && "
        "rm *.tar.* && "
        "cp -R RepeatMasker*/* {params.dir} && "
        "cd {params.dir} && grep -l -r '#!/u1' * | xargs -I '%' sed -i 's|/u1/local/bin/perl|/usr/bin/env perl|' % && "
        "cp RepeatMaskerConfig.tmpl RepeatMaskerConfig.pm && chmod -x RepeatMaskerConfig.pm && "
        'sed -i \'s|DEFAULT_SEARCH_ENGINE\s\+=.*|DEFAULT_SEARCH_ENGINE = "ncbi";|\' RepeatMaskerConfig.pm && '
        'sed -i \'s|TRF_PRGM\s\+=.*|TRF_PRGM = "' + "{input.TRF}" + '";|\' RepeatMaskerConfig.pm && '
        'sed -i \'s|RMBLAST_DIR\s\+=.*|RMBLAST_DIR = "' + "{params.RMBLAST}" + '";|\' RepeatMaskerConfig.pm && '
#        "cp {input.REPBASE}Libraries/RepeatMaskerLib.embl Libraries/ ; "
        "echo 'export PATH='`readlink -f {params.dir}`':$PATH' >> {rules.maker_bashrc.output} && "
        "cd / && rm -rf $TEMPDIR && "
        "touch {output}"

############################################################

fastq-dump --outdir Downloads --split-files ~/ncbi/public/sra/SRR848973.sra
