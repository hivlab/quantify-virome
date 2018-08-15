
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
sample_ids = samples.index.values.tolist()
n_files = config["split_fasta"]["n_files"]

## Target rule
rule all:
    input:
      expand("{sample}/11_bwa_mem/mapped.{n}.bam", sample = sample_ids, n = n_files),
      expand("{sample}/12b_unmapped_masked/RefGenome_unmapped.{n}.masked.fa", sample = sample_ids, n = n_files),
      expand("{sample}/13_megablast/megablast.{n}.xml", sample = sample_ids, n = n_files),
      expand("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.non-viral.out", sample = sample_ids, n = n_files),
      expand("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.unmapped.fa", sample = sample_ids, n = n_files),
      expand("{sample}/15_blast_virusnt/blast_virusnt.{n}.xml", sample = sample_ids, n = n_files),
      expand("{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.known-viral.out", sample = sample_ids, n = n_files),
      expand("{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.unmapped.fa", sample = sample_ids, n = n_files),
      expand("{sample}/reports/taxonomy_report.html", sample = sample_ids)

## Load rules
include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
