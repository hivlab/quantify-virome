
include: "rules/common.smk"

## Target rules
rule all:
    input:
      expand([
      "output/16_blastntvirus_parsed/{sample}_virusnt_blast_{n}_known-viral.out",
      "output/16_blastntvirus_parsed/{sample}_virusnt_blast_{n}_unmapped.fa",
      "output/reports/{sample}_taxonomy_report.html"],
      sample = sample_ids,
      n = list(range(1, n_files + 1, 1)))

## Modules
include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
