
include: "rules/common.smk"

## Target rules
rule all:
    input:
      expand(["output/{sample}_refgenome_filtered_{n}_non-viral.xml",
      "output/{sample}_refgenome_filtered_{n}_unmapped.fa"],
      sample = sample_ids,
      n = list(range(1, n_files + 1, 1)))

## Modules
include: "rules/munge.smk"
include: "rules/cd-hit.smk"
include: "rules/mask.smk"
include: "rules/refgenomefilter.smk"
