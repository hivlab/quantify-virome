
include: "rules/common.smk"

## Target rules
rule all:
    input:
      expand([
      "output/blast/{sample}_megablast_{n}.xml",
      "output/blast/{sample}_blast_virusnt_{n}.xml",
      "output/{sample}_virusnt_blast_{n}_known-viral.out",
      "output/{sample}_virusnt_blast_{n}_unmapped.fa"],
      sample = sample_ids,
      n = list(range(1, n_files + 1, 1)))

## Modules
include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
