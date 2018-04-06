
## Use os.path to update file paths from config file
import os.path
import pandas as pd
shell.executable("bash")

## Load configuration file with sample and path info
configfile: "config.yaml"
# samples = pd.read_csv(config["samples"], index_col = "sample")

## Target rule
rule all:
    input:
      expand(os.path.join(config["outdir"], "repeatmasker_good/{sample}.repeatmasker.goodseq.{mask}.{n}.fa"), sample = config["samples"], n = [1,2,3], mask = ["masked", "unmasked"])

include: "rules/munge.smk"
include: "rules/mask.smk"
