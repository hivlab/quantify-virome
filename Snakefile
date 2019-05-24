
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
from snakemake.remote.zenodo import RemoteProvider as ZENRemoteProvider
from snakemake.utils import validate

# Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")
SAMPLES = pd.read_csv(config["samples"], sep = "\s+").set_index("sample", drop=False)
validate(SAMPLES, "schemas/samples.schema.yaml")
SAMPLE_IDS = SAMPLES.index.values.tolist()
N_FILES = config["split_fasta"]["n_files"]
N = list(range(1, N_FILES + 1, 1))

# Create slurm logs dir
if not os.path.exists("logs/slurm"):
    os.makedirs("logs/slurm")

# Setup Zenodo RemoteProvider
ZEN = ZENRemoteProvider()

# Main output files and target rules
RESULTS = ["phages", "phages_viruses", "non_viral"]
TAXONOMY = expand("taxonomy/{file}.csv",
                file = ["names", "nodes", "division"])
STATS = expand(["stats/{sample}_preprocess.tsv",
                "stats/{sample}_blast.tsv",
                "stats/{sample}_refgenome_stats.txt"],
                sample = SAMPLE_IDS) + expand("stats/{sample}_refbac_stats_{n}.txt",
                sample = SAMPLE_IDS, n = N)
OUTPUTS = expand("results/{sample}_{result}_{n}.csv",
                sample = SAMPLE_IDS, n = N, result = RESULTS) + expand("results/{sample}_unassigned_{n}.fa",
                sample = SAMPLE_IDS, n = N) + TAXONOMY + STATS

# Remote outputs
if config["zenodo"]["deposition_id"]:
    ZENOUTPUTS = [ZEN.remote(expand("{deposition_id}/files/results/{sample}_{result}.csv.tar.gz", deposition_id = config["zenodo"]["deposition_id"], sample = SAMPLE_IDS, result = RESULTS)),
    ZEN.remote(expand("{deposition_id}/files/results/{sample}_unassigned.fa.tar.gz", deposition_id = config["zenodo"]["deposition_id"], sample = SAMPLE_IDS))]

rule all:
    input:
        OUTPUTS, ZENOUTPUTS if config["zenodo"]["deposition_id"] else OUTPUTS

# Modules
include: "rules/preprocess.smk"
include: "rules/blast.smk"
