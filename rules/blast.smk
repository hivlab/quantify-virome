
## MegaBlast against reference genome to remove host sequences [13]
rule megablast_ref_genome:
    input:
      db = config["ref_genome"],
      query = rules.unmapped_masked.output
    output:
      temp("output/13_megablast/{sample}_megablast_{n}.xml")
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
      "output/14_megablast_parsed/{sample}_refgenome_megablast_{n}_non-viral.out",
      "output/14_megablast_parsed/{sample}_refgenome_megablast_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-10
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"

## Blast against virus database [15]
rule blastn_virus_nt:
    input:
      db = config["virus_nt"],
      query = "output/14_megablast_parsed/{sample}_refgenome_megablast_{n}_unmapped.fa"
    output:
      out = "output/15_blast_virusnt/{sample}_blast_virusnt_{n}.xml"
    params:
      task = "blastn",
      show_gis = True,
      evalue = 1e-4,
      db_soft_mask = 100,
      num_threads = 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blastn_virus_db.py"

## Filter blastn records for the cutoff value [16]
rule parse_virusntblast:
    input:
      rules.blastn_virus_nt.output.out,
      "output/14_megablast_parsed/{sample}_refgenome_megablast_{n}_unmapped.fa"
    output:
      "output/16_blastntvirus_parsed/{sample}_virusnt_blast_{n}_known-viral.out",
      "output/16_blastntvirus_parsed/{sample}_virusnt_blast_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"

# Download taxonomy names [17a]
rule download_taxonomy:
    output:
      names = "taxonomy/names.csv",
      nodes = "taxonomy/nodes.csv"
    conda:
      "../envs/tidyverse.yml"
    script:
      "../scripts/download_taxonomy_names.R"

sample_ids, n_ids = glob_wildcards("output/16_blastntvirus_parsed/{sample}_virusnt_blast_{n}_known-viral.out")

# Add taxonomy to virus nt blast [17b]
rule virus_nt_taxonomy:
    input:
      known = expand("output/16_blastntvirus_parsed/{sample}_virusnt_blast_{n}_known-viral.out",
               sample = sample_ids, n = n_ids),
      vhunter = config["vhunter"],
      names = "taxonomy/names.csv",
      nodes = "taxonomy/nodes.csv"
    output:
      "output/17_virus_nt_taxonomy/{sample}_known_taxa.csv"
    conda:
      "../envs/tidyverse.yml"
    script:
      "../scripts/munge_taxonomy.R"

# Taxonomy report to virus nt blast [17c]
rule virus_nt_taxonomy_report:
    input:
      rules.virus_nt_taxonomy.output,
      names = "taxonomy/names.csv"
    output:
      "output/reports/{sample}_taxonomy_report.html"
    params:
      lambda wildcards: wildcards.sample
    conda:
      "../envs/tidyverse.yml"
    script:
      "../scripts/taxonomy_report.Rmd"
