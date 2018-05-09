
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

## Use os.path to update file paths from config file
import os.path
import pandas as pd
shell.executable("bash")

## Load configuration file with sample and path info
configfile: "config.yml"
samples = pd.read_table(config["samples"], sep = "\s+", index_col = "sample", dtype = str)

## Target rule
rule all:
    input:
      expand(os.path.join(config["outdir"], "{sample}/14_megablast_parsed/megablast_parsed.{n}.{ext}"), sample = "I1164_12629_Harvard_SIV_196_06_2_24_12_mini", n = [1, 2], ext = ["out"])

include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
