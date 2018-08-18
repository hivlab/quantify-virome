
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
      expand([os.path.join(config["outdir"], "{sample}/01_fastp/pair1.truncated.gz"),
      os.path.join(config["outdir"], "{sample}/01_fastp/pair1.truncated.gz"),
      os.path.join(config["outdir"], "{sample}/02_stitched/join.fq.gz"),
      os.path.join(config["outdir"], "{sample}/02_stitched/un1.fq.gz"),
      os.path.join(config["outdir"], "{sample}/02_stitched/un2.fq.gz"),
      os.path.join(config["outdir"], "{sample}/03_merged/stitched.merged.fq.gz"),
      os.path.join(config["outdir"], "{sample}/04_fasta/stitched.merged.fasta"),
      os.path.join(config["outdir"], "{sample}/05_cdhit/merged.cdhit.fa"),
      os.path.join(config["outdir"], "{sample}/06_tantan/cdhit.tantan.fa"),
      os.path.join(config["outdir"], "{sample}/07_tantan_good/tantan.goodseq.fa"),
      os.path.join(config["outdir"], "{sample}/17_virus_nt_taxonomy/known_taxa.csv"),
      os.path.join(config["outdir"], "{sample}/17_virus_nt_taxonomy/taxonomy_report.html")],
      sample = samples.index.values.tolist()),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/08_split_fasta/tantan.goodseq.{{n}}.fa"),
      sample = samples.index.values.tolist()))

include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
