
include: "rules/common.smk"

## Target rules
rule all:
    input:
      "{sample}/11_bwa_mem/mapped.{n}.bam",
      "{sample}/12b_unmapped_masked/RefGenome_unmapped.{n}.masked.fa",
      "{sample}/13_megablast/megablast.{n}.xml",
      "{sample}/14_megablast_parsed/RefGenome_megablast.{n}.non-viral.out",
      "{sample}/14_megablast_parsed/RefGenome_megablast.{n}.unmapped.fa",
      "{sample}/15_blast_virusnt/blast_virusnt.{n}.xml",
      "{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.known-viral.out",
      "{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.unmapped.fa"

## Modules
include: "rules/munge.smk"
include: "rules/mask.smk"
include: "rules/align.smk"
include: "rules/blast.smk"
