
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

include: "rules/common.smk"

## Target rules
rule all:
    input:
      expand([
      "output/blast/{sample}_megablast_{n}.xml",
      "output/blast/{sample}_blast_virusnt_{n}.xml"],
      sample = sample_ids,
      n = list(range(1, n_files + 1, 1)))

## Modules
include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
