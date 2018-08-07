
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
      dynamic(expand("{sample}/11_bwa_mem/mapped.{{n}}.bam", sample = samples.index.values.tolist())),
      expand("{sample}/12b_unmapped_masked/RefGenome_unmapped.{n}.masked.fa", sample = samples.index.values.tolist()),
      expand("{sample}/13_megablast/megablast.{n}.xml", sample = samples.index.values.tolist()),
      expand("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.non-viral.out", sample = samples.index.values.tolist()),
      expand("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.unmapped.fa", sample = samples.index.values.tolist()),
      expand("{sample}/15_blast_virusnt/blast_virusnt.{n}.xml", sample = samples.index.values.tolist()),
      expand("{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.known-viral.out", sample = samples.index.values.tolist()),
      expand("{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.unmapped.fa", sample = samples.index.values.tolist()),
      expand("reports/{sample}/taxonomy_report.html", sample = samples.index.values.tolist())

## Setup report

## report: "report/workflow.rst"

## Load rules
include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
