
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

include: "rules/common.smk"

## Target rules
rule all:
    input:
      expand([
      "results/{sample}_phages_{n}.csv",
      "results/{sample}_unassigned_blasted_{n}.fa",
      "results/{sample}_phages_blasted_{n}.csv",
      "results/{sample}_viruses_blasted_{n}.csv"
      ],
      sample = sample_ids,
      n = list(range(1, n_files + 1, 1))),
      expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])

## Modules
include: "rules/munge.smk"
include: "rules/cd-hit.smk"
include: "rules/mask.smk"
include: "rules/refgenomefilter.smk"
include: "rules/blast.smk"
