__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

## Load libraries
import os
import json
import glob
import pandas as pd
from snakemake.remote.zenodo import RemoteProvider as ZENRemoteProvider
from snakemake.utils import validate

ZEN = ZENRemoteProvider()

## Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")
SAMPLES = pd.read_table(config["samples"], sep = "\s+").set_index("sample", drop=False)
validate(SAMPLES, "schemas/samples.schema.yaml")
SAMPLE_IDS = SAMPLES.index.values.tolist()
N_FILES = config["split_fasta"]["n_files"]
N = list(range(1, N_FILES + 1, 1))

localrules: all, upload_stats

rule all:
   input: ZEN.remote(expand("1488086/files/{sample}_stats.json", sample = SAMPLE_IDS))

rule collect_stats:
   input:
      fastp = "munge/{sample}_fastp_report.json",
      cdhit = "logs/{sample}_cdhit.log",
      rm = expand("mask/{{sample}}_repeatmasker_{n}.fa.tbl", n = N),
      rm_good = expand("mask/{{sample}}_unmaskedgood_{n}.fa", n = N),
      refgenome_unmapped = expand("refgenomefilter/{{sample}}_refgenome_unmapped_{n}.fa", n = N),
      refgenome_filtered = expand("refgenomefilter/{{sample}}_refgenome_filtered_{n}_known-host.tsv", n = N)
   output: 
      outfile = "results/{sample}_stats.json" 
   script:
      "../scripts/collect_stats.py"

if config["zenodo"]["deposition_id"]:
   rule upload_stats:
      input: rules.collect_stats.output
      output: ZEN.remote("1488086/files/{sample}_stats.json")
      shell: "cp {input} {output}"
