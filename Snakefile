
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

subworkflow virome:
    snakefile: "virome.snakefile"

include: "rules/common.smk"
include: "rules/taxonomy.smk"

rule all:
  input: "taxonomy/names.csv", "taxonomy/nodes.csv", expand(["output/reports/{sample}_known_taxa.csv", "output/reports/{sample}_taxonomy_report.html"], sample = sample_ids)
