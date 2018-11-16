
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

include: "rules/common.smk"

## Main output files
file_ids = list(range(1, n_files + 1, 1))
outputs = expand(["results/{sample}_phages_{n}.csv", "results/{sample}_unassigned_{n}.fa", "results/{sample}_phages_blasted_{n}.csv", "results/{sample}_viruses_blasted_{n}.csv"], sample = sample_ids, n = file_ids) + expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])

## Target rules
rule all:
    input:
        outputs, expand("results/{sample}_phages.csv.tar.gz", sample = sample_ids) if config["zenodo"]["deposition_id"] else outputs

## Modules
include: "rules/munge.smk"
include: "rules/cd-hit.smk"
include: "rules/mask.smk"
include: "rules/refgenomefilter.smk"
include: "rules/blast.smk"
