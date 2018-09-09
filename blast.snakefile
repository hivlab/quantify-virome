
subworkflow preprocessing:
    snakefile: "preprocessing.snakefile"

include: "rules/common.smk"

rule all:
    input:
      expand([
      "output/blast/{sample}_blastn_virus_{n}.xml",
      "output/{sample}_blastn_virus_{n}_known-viral.xml",
      "output/{sample}_blastn_virus_{n}_known-viral.fa",
      "output/{sample}_blastn_virus_{n}_unmapped.fa",
      "output/blast/{sample}_blastx_virus_{n}.xml",
      "output/{sample}_blastx_virus_{n}_known-viral.xml",
      "output/{sample}_blastx_virus_{n}_known-viral.fa",
      "output/{sample}_blastx_virus_{n}_unmapped.fa",
      "output/{sample}_known-viral_{n}_unmasked.fa"
      ],
      sample = sample_ids,
      n = list(range(1, n_files + 1, 1)))

include: "rules/blast.smk"
