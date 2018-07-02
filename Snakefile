
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

## Use os.path to update file paths from config file
import os.path
from os import listdir
import glob
import pandas as pd
shell.executable("bash")

## Load configuration file with sample and path info
configfile: "config.yml"
samples = pd.read_table(config["samples"], sep = "\s+", index_col = "sample", dtype = str)

## Target rule
rule all:
    input:
      dynamic(expand(os.path.join(config["outdir"], "{sample}/10_repeatmasker_good/masked.{{n}}.fa"), sample = samples.index.values.tolist())), dynamic(expand(os.path.join(config["outdir"], "{sample}/10_repeatmasker_good/unmasked.{{n}}.fa"), sample = samples.index.values.tolist()))

include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
