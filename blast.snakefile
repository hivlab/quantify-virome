
subworkflow preprocessing:
    snakefile: "preprocessing.snakefile"

include: "rules/common.smk"

rule all:
    input:
      expand([
      "output/{sample}_phages_{n}.csv",
      "output/{sample}_blastn_nt_{n}_mapped.xml",
      "output/{sample}_blastn_nt_{n}_unmapped.fa",
      "output/{sample}_blastx_nr_{n}_mapped.xml",
      "output/{sample}_blastx_nr_{n}_unmapped.fa"
      ], sample = sample_ids, n = list(range(1, n_files + 1, 1))),
      "taxonomy/names.csv", "taxonomy/nodes.csv", "taxonomy/division.csv"

include: "rules/blast.smk"
