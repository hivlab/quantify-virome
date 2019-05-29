
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
shell.executable("bash")

# Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")
RUNS = pd.read_csv(config["samples"], sep = "\s+")
index = pd.MultiIndex.from_frame(RUNS[["group", "run"]])
RUNS = RUNS.set_index(index)
validate(RUNS, "schemas/samples.schema.yaml")
GROUP_IDS,RUN_IDS = list(zip(*RUNS.index.to_list()))
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
STATS = expand(["stats/{run}_preprocess.tsv",
                "stats/{run}_blast.tsv",
                "stats/{run}_refgenome_stats.txt"],
                run = RUN_IDS) + expand("stats/{run}_refbac_stats_{n}.txt",
                run = RUN_IDS, n = N)
OUTPUTS = expand("results/{run}_{result}_{n}.csv",
                run = RUN_IDS, n = N, result = RESULTS) + expand("results/{run}_unassigned_{n}.fa",
                run = RUN_IDS, n = N) + TAXONOMY + STATS

# Remote outputs
if config["zenodo"]["deposition_id"]:
    ZENOUTPUTS = [ZEN.remote(expand("{deposition_id}/files/results/{run}_{result}.csv.tar.gz", deposition_id = config["zenodo"]["deposition_id"], run = RUN_IDS, result = RESULTS)),
    ZEN.remote(expand("{deposition_id}/files/results/{run}_unassigned.fa.tar.gz", deposition_id = config["zenodo"]["deposition_id"], run = RUN_IDS))]

rule all:
    input:
        OUTPUTS, ZENOUTPUTS if config["zenodo"]["deposition_id"] else OUTPUTS

# Modules
include: "rules/preprocess.smk"
include: "rules/blast.smk"
