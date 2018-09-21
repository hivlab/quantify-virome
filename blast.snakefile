
subworkflow preprocessing:
    snakefile: "preprocessing.snakefile"

include: "rules/common.smk"

rule all:
    input:
      expand([
      "output/{sample}_phages_{n}.csv",
      "output/{sample}_phages_blasted_{n}.csv",
      "output/{sample}_viruses_blasted_{n}.csv"
      ], sample = sample_ids, n = list(range(1, n_files + 1, 1))),
      "taxonomy/names.csv", "taxonomy/nodes.csv", "taxonomy/division.csv"

include: "rules/blast.smk"
