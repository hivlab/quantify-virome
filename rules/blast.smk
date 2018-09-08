

## Blast against NT virus database
rule blastn_virus:
    input:
      db = config["virus_nt"],
      query = preprocessing("output/{sample}_refgenome_filtered_{n}_unmapped.fa")
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

## Filter blastn records for the cutoff value
rule parse_blastn_virus:
    input:
      rules.blastn_virus.output,
      rules.blastn_virus.input.query
    output:
      known = "output/{sample}_blastn_virus_{n}_known-viral.xml",
      unmapped = "output/{sample}_blastn_virus_{n}_unmapped.fa"
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
      query = rules.parse_blastn_virus.output.unmapped
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
      rules.blastx_virus.input.query
    output:
      known = "output/{sample}_blastx_virus_{n}_known-viral.xml",
      unmapped = "output/{sample}_blastx_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-3
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"
