
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
      expand(os.path.join(config["outdir"], "{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.{ext}"), sample = "I1164_12629_Harvard_SIV_196_06_2_24_12_mini", n = [1, 2], ext = ["bam"]),
      expand(os.path.join(config["outdir"], "{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.{ext}"), sample = "test_seq_001", n = 1, ext = ["bam"])

include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
