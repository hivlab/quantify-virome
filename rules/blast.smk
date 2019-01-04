
## Prepare names and nodes tables for taxonomy annotation
rule prepare_taxonomy_data:
  input: config["names"], config["nodes"], config["division"]
  output:
      expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
  conda:
    "../envs/tidyverse.yaml"
  script:
    "../scripts/prepare_taxonomy_data.R"

## Blast input, output, and params keys must match commandline blast option names https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a

## Blast against NT virus database
rule blastn_virus:
    input:
      query = rules.parse_megablast.output.unmapped,
      db = config["virus_nt"]
    output:
      out = "blast/{sample}_blastn_virus_{n,\d+}.tsv"
    params:
      task = "blastn",
      evalue = config["blastn_virus"]["evalue"],
      db_soft_mask = config["blastn_virus"]["db_soft_mask"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

## Filter blastn records for the cutoff value
rule parse_blastn_virus:
    input:
      query = rules.parse_megablast.output.unmapped,
      blast_result = rules.blastn_virus.output.out
    output:
      mapped = "blast/{sample}_blastn_virus_{n,\d+}_known-viral.tsv",
      unmapped = "blast/{sample}_blastn_virus_{n,\d+}_unmapped.fa"
    params:
      e_cutoff = 1e-5,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

## Blastx unmapped sequences against NR virus database
rule blastx_virus:
    input:
      query = rules.parse_blastn_virus.output.unmapped,
      db = config["virus_nr"]
    output:
      out = "blast/{sample}_blastx_virus_{n}.tsv"
    params:
      word_size = 6,
      evalue = config["blastx_virus"]["evalue"],
      db_soft_mask = config["blastx_virus"]["db_soft_mask"],
      max_hsps = config["blastx_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

## Filter blastn records for the cutoff value
rule parse_blastx_virus:
    input:
      query = rules.blastx_virus.input.query,
      blast_result = rules.blastx_virus.output.out
    output:
      mapped = "blast/{sample}_blastx_virus_{n}_known-viral.tsv",
      unmapped = "blast/{sample}_blastx_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-3,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

## Filter out phage sequences
rule filter_viruses:
  input:
    [rules.parse_blastn_virus.output.mapped,
    rules.parse_blastx_virus.output.mapped] if config["run_blastx"] else rules.parse_blastn_virus.output.mapped,
    taxdb = config["vhunter"],
    nodes = "taxonomy/nodes.csv"
  output:
    phages = "results/{sample}_phages_{n}.csv",
    viruses = "blast/{sample}_candidate_viruses_{n}.csv"
  conda:
    "../envs/tidyverse.yaml"
  script:
    "../scripts/filter_viruses.R"

## Get unmasked candidate viral sequences
rule unmasked_viral:
    input:
      rules.filter_viruses.output.viruses,
      rules.refgenome_unmapped.output.fa
    output:
      "blast/{sample}_candidate_viruses_{n}_unmasked.fa"
    conda:
      "../envs/biopython.yaml"
    script:
      "../scripts/unmasked_viral.py"

## Map against bacterial genomes
rule bwa_map_refbac:
    input:
      config["ref_bacteria"],
      [rules.unmasked_viral.output]
    output:
      "blast/{sample}_bac_mapped_{n}.sam"
    log:
      "logs/{sample}_bwa_map_refbac_{n}.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yaml"
    shell:
      "bwa mem -k 15 -t {threads} {input} > {output} 2> {log}"

## Extract unmapped reads
rule refbac_unmapped:
    input: rules.bwa_map_refbac.output
    output:
      bam = "blast/{sample}_bac_unmapped_{n}.bam",
      fq = "blast/{sample}_bac_unmapped_{n}.fq",
      fa = "blast/{sample}_bac_unmapped_{n}.fa"
    conda:
      "../envs/bwa-sam-bed.yaml"
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
      rules.repeatmasker_good.output.masked_filt
    output:
      "blast/{sample}_bac_unmapped_{n}_masked.fa"
    conda:
      "../envs/biopython.yaml"
    script:
      "../scripts/unmapped_masked_ids.py"

## MegaBlast against NT to remove host sequences
rule megablast_nt:
    input:
      query = rules.refbac_unmapped_masked.output,
      db = config["nt"]
    output:
      out = "blast/{sample}_megablast_nt_{n,\d+}.tsv"
    params:
      task = "megablast",
      evalue = config["megablast_nt"]["evalue"],
      word_size = config["megablast_nt"]["word_size"],
      max_hsps = config["megablast_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

## Filter megablast records for the cutoff value
rule parse_megablast_nt:
    input:
      query = rules.refbac_unmapped_masked.output,
      blast_result = rules.megablast_nt.output.out
    output:
      mapped = "blast/{sample}_megablast_nt_{n,\d+}_mapped.tsv",
      unmapped = "blast/{sample}_megablast_nt_{n,\d+}_unmapped.fa"
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

## Blastn against NT database
rule blastn_nt:
    input:
      query = rules.parse_megablast_nt.output.unmapped,
      db = config["nt"]
    output:
      out = "blast/{sample}_blastn_nt_{n,\d+}.tsv"
    params:
      task = "blastn",
      evalue = config["blastn_nt"]["evalue"],
      max_hsps = config["blastn_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

## Filter blastn records for the cutoff value
rule parse_blastn_nt:
    input:
      query = rules.blastn_nt.input.query,
      blast_result = rules.blastn_nt.output.out
    output:
      mapped = "blast/{sample}_blastn_nt_{n,\d+}_mapped.tsv",
      unmapped = "blast/{sample}_blastn_nt_{n,\d+}_unmapped.fa" if config["run_blastx"] else "results/{sample}_unassigned_{n,\d+}.fa"
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

## Blastx unmapped sequences against NR virus database
rule blastx_nr:
    input:
      query = rules.parse_blastn_nt.output.unmapped,
      db = config["nr"]
    output:
      out = "blast/{sample}_blastx_nr_{n,\d+}.tsv"
    params:
      word_size = 6,
      evalue = config["blastx_nr"]["evalue"],
      max_hsps = config["blastx_nr"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

## Filter blastn records for the cutoff value
rule parse_blastx_nr:
    input:
      query = rules.blastx_nr.input.query,
      blast_result = rules.blastx_nr.output.out
    output:
      mapped = "blast/{sample}_blastx_nr_{n,\d+}_mapped.tsv",
      unmapped = "results/{sample}_unassigned_{n,\d+}.fa" if config["run_blastx"] else "{sample}_None_{n}"
    params:
      e_cutoff = 1e-3,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

## Filter out virus and phage sequences
rule filter_blasted_viruses:
  input:
    [rules.parse_blastn_nt.output.mapped, rules.parse_blastx_nr.output.mapped] if config["run_blastx"] else rules.parse_blastn_nt.output.mapped,
    taxdb = config["vhunter"],
    nodes = "taxonomy/nodes.csv"
  output:
    phages = "results/{sample}_phages_blasted_{n}.csv",
    viruses = "results/{sample}_viruses_blasted_{n}.csv"
  conda:
    "../envs/tidyverse.yaml"
  script:
    "../scripts/filter_viruses.R"

## Upload results to Zenodo
if config["zenodo"]["deposition_id"]:
    rule upload:
        input:
          expand("results/{{sample}}_{{result}}_{n}.{{ext}}", n = N)
        output:
          "results/{sample}_{result}.{ext}.tar.gz"
        params:
          config["zenodo"]["deposition_id"]
        wrapper:
          config["wrappers"]["upload"]
