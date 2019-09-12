
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import json
import glob
import pandas as pd
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.utils import validate
shell.executable("bash")

# Load configuration file with sample and path info
configfile: "config.yaml"
#validate(config, "schemas/config.schema.yaml")

# Load runs and groups
RUNS = pd.read_csv(config["samples"], sep = "\s+").set_index("run", drop = False)
validate(RUNS, "schemas/samples.schema.yaml")
RUN_IDS = RUNS.index.tolist()
N_FILES = config["split_fasta"]["n_files"]
N = list(range(1, N_FILES + 1, 1))

# Create slurm logs dir
if not os.path.exists("logs/slurm"):
    os.makedirs("logs/slurm")

wildcard_constraints:
    run = "[a-zA-Z0-9]+",
    n = "\d+"

# Main output files
RESULTS = ["phages.csv", "phages-viruses.csv", "non-viral.csv", "unassigned.fa"]
BLASTV = ["blastn-virus", "blastx-virus"] if config["run_blastx"] else ["blastn-virus"]
BLASTNR = ["megablast-nt", "blastn-nt", "blastx-nr"] if config["run_blastx"] else ["megablast-nt", "blastn-nt"]
BLAST = BLASTV + BLASTNR
STATS = expand(["stats/{run}_refgenome-stats.txt", "stats/{run}_preprocess.tsv", "stats/{run}_blast.tsv"], run = RUN_IDS) + expand("stats/{run}_refbac-stats_{n}.txt", run = RUN_IDS, n = N)
OUTPUTS = expand("results/{run}_{result}", run = RUN_IDS, result = RESULTS) + STATS

# Remote outputs
if config["zenodo"]["deposition_id"]:
    # Load zenodo remote provider module
    from snakemake.remote.zenodo import RemoteProvider as ZENRemoteProvider
    # Setup Zenodo RemoteProvider
    ZEN = ZENRemoteProvider(deposition = config["zenodo"]["deposition_id"], access_token = os.environ["ZENODO_PAT"])
    # Append uploads
    ZENOUTPUTS = ZEN.remote(expand(["results/{run}_counts.tgz", "stats/{run}_stats.tgz"], run = RUN_IDS))
    OUTPUTS = OUTPUTS + ZENOUTPUTS

rule all:
    input:
        OUTPUTS

# Path to reference genomes
REF_GENOME = os.getenv("REF_GENOME_HUMAN")
REF_BACTERIA = os.getenv("REF_BACTERIA")

# Wrappers
BLAST = "https://raw.githubusercontent.com/avilab/virome-wrappers/blast5/blast/query"
PARSE_BLAST = "https://raw.githubusercontent.com/avilab/virome-wrappers/master/blast/parse"
STATS = "https://bitbucket.org/tpall/snakemake-wrappers/raw/e7699c0ae37a999909fb764c91723d46ded7461c/bio/seqkit/stats"

# Rules
include: "rules/preprocess.smk"
include: "rules/blast.smk"
