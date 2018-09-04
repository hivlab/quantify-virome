
localrules: download_taxonomy

# Download taxonomy names [17a]
rule download_taxonomy:
    input:
      FTP.remote("ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz")
    output:
      names_dmp = temp("names.dmp"),
      nodes_dmp = temp("nodes.dmp"),
      names = "taxonomy/names.csv",
      nodes = "taxonomy/nodes.csv"
    params:
      taxdir = "taxonomy"
    conda:
      "../envs/tidyverse.yml"
    script:
      "../scripts/download_taxonomy_names.R"

def get_knownviral(wildcards):
  path = expand("output/{sample}_virusnt_blast_*_known-viral.out", sample = wildcards.sample)
  return glob.glob(*path)

# Add taxonomy to virus nt blast [17b]
rule virus_nt_taxonomy:
    input:
      known = get_knownviral,
      vhunter = config["vhunter"],
      names = rules.download_taxonomy.output.names,
      nodes = rules.download_taxonomy.output.nodes
    output:
      "output/reports/{sample}_known_taxa.csv"
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
