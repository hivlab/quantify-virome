
## MegaBlast against reference genome to remove host sequences [13]
rule megablast_ref_genome:
    input:
      db = config["ref_genome"],
      query = rules.unmapped_masked.output
    output:
      "output/blast/{sample}_megablast_{n}.xml"
    params:
      perc_ident = config["megablast_ref_genome"]["perc_identity"],
      evalue = config["megablast_ref_genome"]["evalue"],
      word_size = config["megablast_ref_genome"]["word_size"],
      num_desc = config["megablast_ref_genome"]["num_descriptions"],
      num_align = config["megablast_ref_genome"]["num_alignments"]
    threads: 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/megablast_ref_genome.py"

## Filter megablast records for the cutoff value [14]
rule parse_megablast:
    input:
      rules.megablast_ref_genome.output,
      rules.unmapped_masked.output
    output:
      "output/{sample}_refgenome_megablast_{n}_non-viral.out",
      "output/{sample}_refgenome_megablast_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-10
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"

## Blast against NT virus database [15]
rule blastn_virus:
    input:
      db = config["virus_nt"],
      query = "output/{sample}_refgenome_megablast_{n}_unmapped.fa"
    output:
      "output/blast/{sample}_blastn_virus_{n}.xml"
    params:
      show_gis = True,
      evalue = 1e-4,
      db_soft_mask = 100
    threads: 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blastn_virus.py"

## Filter blastn records for the cutoff value [16]
rule parse_blastn_virus:
    input:
      rules.blastn_virus.output,
      "output/{sample}_refgenome_megablast_{n}_unmapped.fa"
    output:
      "output/{sample}_blastn_virus_{n}_known-viral.out",
      "output/{sample}_blastn_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"

## Blastx unmapped sequences against NR virus database
rule blastx_virus:
    input:
      db = config["virus_nr"],
      query = "output/{sample}_blastn_virus_{n}_unmapped.fa"
    output:
      "output/blast/{sample}_blastx_virus_{n}.xml"
    params:
      show_gis = True,
      evalue = 1e-2,
      db_soft_mask = 100
    threads: 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blastx_virus.py"

## Filter blastn records for the cutoff value
rule parse_blastx_virus:
    input:
      rules.blastx_virus.output,
      "output/blast/{sample}_blastx_virus_{n}.xml"
    output:
      "output/{sample}_blastx_virus_{n}_known-viral.out",
      "output/{sample}_blastx_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-3
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"
