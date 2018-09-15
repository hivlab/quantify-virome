
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
      "output/{sample}_known-viral_{n}_unmasked.fa",
      "output/{sample}_refbacteria_unmapped_{n}_masked.fa",
      "output/blast/{sample}_nt_filtered_{n}_mapped.xml",
      "output/blast/{sample}_nt_filtered_{n}_unmapped.fa",
      "output/{sample}_blastn_nt_{n}_mapped.xml",
      "output/{sample}_blastn_nt_{n}_unmapped.fa",
      "output/{sample}_blastx_nr_{n}_mapped.xml",
      "output/{sample}_blastx_nr_{n}_unmapped.fa",
      "output/{sample}_phages_{n}.csv",
      "output/{sample}_candidate_viruses_{n}.csv"
      ], sample = sample_ids, n = list(range(1, n_files + 1, 1))),
      "taxonomy/names.csv", "taxonomy/nodes.csv", "taxonomy/division.csv"

include: "rules/blast.smk"
