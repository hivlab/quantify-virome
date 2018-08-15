
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
      os.path.join(config["outdir"], "{sample}/07_tantan_good/tantan.goodseq.fa")],
      sample = samples.index.values.tolist()),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/08_split_fasta/tantan.goodseq.{{n}}.fa"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/09_repeatmasker/tantan.goodseq.{{n}}.fa.masked"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/10_repeatmasker_good/masked.{{n}}.fa"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/10_repeatmasker_good/unmasked.{{n}}.fa"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/11_bwa_mem/mapped.{{n}}.bam"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/12a_unmapped_reads/RefGenome_unmapped.{{n}}.bam"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/12a_unmapped_reads/RefGenome_unmapped.{{n}}.fq"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/12a_unmapped_reads/RefGenome_unmapped.{{n}}.fa"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/12b_unmapped_masked/RefGenome_unmapped.{{n}}.masked.fa"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/13_megablast/megablast.{{n}}.xml"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/14_megablast_parsed/RefGenome_megablast.{{n}}.non-viral.out"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/14_megablast_parsed/RefGenome_megablast.{{n}}.unmapped.fa"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/15_blast_virusnt/blast_virusnt.{{n}}.xml"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/16_blastntvirus_parsed/blastnt_virus.{{n}}.known-viral.out"),
      sample = samples.index.values.tolist())),
      dynamic(expand(os.path.join(config["outdir"], "{sample}/16_blastntvirus_parsed/blastnt_virus.{{n}}.unmapped.fa"),
      sample = samples.index.values.tolist())),
      os.path.join(config["datadir"], "names.csv"),
      os.path.join(config["datadir"], "nodes.csv")

include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
