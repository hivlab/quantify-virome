

def get_knownviral(wildcards):
  path = expand([
      "{sample}_phages_*.csv",
      "{sample}_phages_blasted_*.csv",
      "{sample}_viruses_blasted_*.csv"
      ], sample = wildcards.sample)
  csv = [item for sublist in list(map(glob.glob, path)) for item in sublist]
  return csv

# Taxonomy report to virus nt blast
rule virus_nt_taxonomy_report:
    input:
      get_knownviral,
      "taxonomy/names.csv",
      "taxonomy/division.csv"
    output:
      reports/{sample}_taxonomy_report.html"
    params:
      lambda wildcards: wildcards.sample
    conda:
      "../envs/tidyverse.yml"
    script:
      "../scripts/taxonomy_report.Rmd"
