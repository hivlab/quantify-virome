
## Prepare names and nodes tables for taxonomy annotation
rule prepare_taxonomy_data:
  input: config["names"], config["nodes"], config["division"]
  output:
      expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
  conda:
    "../envs/tidyverse.yml"
  script:
    "../scripts/prepare_taxonomy_data.R"

## Blast input, output, and params keys must match commandline blast option names https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a

## Blast against NT virus database
rule blastn_virus:
    input:
      query = "{sample}_refgenome_filtered_{n}_unmapped.fa"
    output:
      out = protected("blast/{sample}_blastn_virus_{n}.xml")
    params:
      db = config["virus_nt"],
      task = "blastn",
      evalue = config["blastn_virus"]["evalue"],
      db_soft_mask = config["blastn_virus"]["db_soft_mask"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastn_virus:
    input:
      rules.blastn_virus.output.out,
      "{sample}_refgenome_filtered_{n}_unmapped.fa"
    output:
      known_xml = temp("{sample}_blastn_virus_{n}_known-viral.xml"),
      unmapped = temp("{sample}_blastn_virus_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Blastx unmapped sequences against NR virus database
rule blastx_virus:
    input:
      query = rules.parse_blastn_virus.output.unmapped
    output:
      out = protected("blast/{sample}_blastx_virus_{n}.xml")
    params:
      db = config["virus_nr"],
      evalue = config["blastx_virus"]["evalue"],
      db_soft_mask = config["blastx_virus"]["db_soft_mask"],
      max_hsps = config["blastx_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastx_virus:
    input:
      rules.blastx_virus.output.out,
      rules.blastx_virus.input.query
    output:
      known_xml = temp("{sample}_blastx_virus_{n}_known-viral.xml"),
      unmapped = temp("{sample}_blastx_virus_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-3
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Filter out phage sequences
rule filter_viruses:
  input:
    rules.parse_blastn_virus.output.known_xml,
    rules.parse_blastx_virus.output.known_xml,
    taxdb = config["vhunter"],
    nodes = "taxonomy/nodes.csv"
  output:
    phages = protected("{sample}_phages_{n}.csv"),
    viruses = protected("{sample}_candidate_viruses_{n}.csv")
  conda:
    "../envs/tidyverse.yml"
  script:
    "../scripts/filter_viruses.R"

## Get unmasked candidate viral sequences
rule unmasked_viral:
    input:
      rules.filter_viruses.output.viruses,
      "{sample}_refgenome_unmapped_{n}.fa"
    output:
      "{sample}_candidate_viruses_{n}_unmasked.fa"
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/unmasked_viral.py"

## Map against bacterial genomes
rule map_refbac:
    input:
      config["ref_bacteria"],
      ["{sample}_candidate_viruses_{n}_unmasked.fa"]
    output:
      temp("{sample}_bac_mapped_{n}.sam")
    log:
      "logs/{sample}_map_refbac_{n}.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      "bwa mem -k 15 -t {threads} {input} > {output} 2> {log}"

## Extract unmapped reads
rule refbac_unmapped:
    input: rules.map_refbac.output
    output:
      bam = temp("{sample}_bac_unmapped_{n}.bam"),
      fq = temp("{sample}_bac_unmapped_{n}.fq"),
      fa = "{sample}_bac_unmapped_{n}.fa"
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
      samtools view -b -S -f 4 {input} > {output.bam}
      bedtools bamtofastq -i {output.bam} -fq {output.fq}
      cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

## Subset repeatmasker masked reads using unmapped ids
rule refbac_unmapped_masked:
    input:
      rules.refbac_unmapped.output.fa,
      "{sample}_repmaskedgood_{n}.fa"
    output:
      temp("{sample}_bac_unmapped_{n}_masked.fa")
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/unmapped_masked_ids.py"

## MegaBlast against NT to remove host sequences
rule megablast_nt:
    input:
      query = rules.refbac_unmapped_masked.output
    output:
      out = "blast/{sample}_megablast_nt_{n}.xml"
    params:
      db = config["nt"],
      task = "megablast",
      evalue = config["megablast_nt"]["evalue"],
      word_size = config["megablast_nt"]["word_size"],
      max_hsps = config["megablast_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter megablast records for the cutoff value
rule parse_megablast_nt:
    input:
      rules.megablast_nt.output.out,
      rules.refbac_unmapped_masked.output
    output:
      known_xml = temp("{sample}_nt_filtered_{n}_mapped.xml"),
      unmapped = temp("{sample}_nt_filtered_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Blastn against NT database
rule blastn_nt:
    input:
      query = rules.parse_megablast_nt.output.unmapped
    output:
      out = "blast/{sample}_blastn_nt_{n}.xml"
    params:
      db = config["nt"],
      task = "blastn",
      evalue = config["blastn_nt"]["evalue"],
      max_hsps = config["blastn_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastn_nt:
    input:
      rules.blastn_nt.output.out,
      rules.blastn_nt.input.query
    output:
      known_xml = temp("{sample}_blastn_nt_{n}_mapped.xml"),
      unmapped = temp("{sample}_blastn_nt_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Blastx unmapped sequences against NR virus database
rule blastx_nr:
    input:
      query = rules.parse_blastn_nt.output.unmapped
    output:
      out = "blast/{sample}_blastx_nr_{n}.xml"
    params:
      db = config["nr"],
      evalue = config["blastx_nr"]["evalue"],
      max_hsps = config["blastx_nr"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastx_nr:
    input:
      rules.blastx_nr.output.out,
      rules.blastx_nr.input.query
    output:
      known_xml = temp("{sample}_blastx_nr_{n}_mapped.xml"),
      unassigned = temp("{sample}_blastx_nr_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-3
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Filter out phage sequences
rule filter_blasted_viruses:
  input:
    "{sample}_blastn_nt_{n}_mapped.xml",
    "{sample}_blastx_nr_{n}_mapped.xml",
    taxdb = config["vhunter"],
    nodes = "taxonomy/nodes.csv"
  output:
    phages = "{sample}_phages_blasted_{n}.csv",
    viruses = "{sample}_viruses_blasted_{n}.csv"
  conda:
    "../envs/tidyverse.yml"
  script:
    "../scripts/filter_viruses.R"
