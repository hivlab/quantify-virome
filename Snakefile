
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

## Use os.path to update file paths from config file
import os.path
from pathlib import Path
import pandas as pd
shell.executable("bash")

## Load configuration file with sample and path info
configfile: "config.yml"
samples = pd.read_table(config["samples"], sep = "\s+", index_col = "sample", dtype = str)

## Target rule
rule all:
    input:
      expand(os.path.join(config["outdir"], "{sample}/17_virus_nt_taxonomy/known_taxa.csv"), sample = "I1164_12629_Harvard_SIV_196_06_2_24_12"),
      expand(os.path.join(config["outdir"], "{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.{ext}"), sample = "I1164_12629_Harvard_SIV_196_06_2_24_12", n = range(1, 24), ext = ["known-viral.out", "unmapped.fa"])


include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
