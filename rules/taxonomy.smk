

def get_knownviral(wildcards):
  path = expand("output/{sample}_virusnt_blast_*_known-viral.out", sample = wildcards.sample)
  return glob.glob(*path)

# Add taxonomy to virus nt blast
rule virus_nt_taxonomy:
    input:
      known = get_knownviral,
      vhunter = config["vhunter"],
      names = config["names"],
      nodes = config["nodes"]
    output:
      "output/reports/{sample}_known_taxa.csv"
    conda:
      "../envs/tidyverse.yml"
    script:
      "../scripts/munge_taxonomy.R"

# Taxonomy report to virus nt blast [35+]
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
