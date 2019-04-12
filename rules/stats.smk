
rule collect_stats:
   input:
      fastp = "munge/{sample}_fastp_report.json",
      cdhit = "logs/{sample}_cdhit.log",
      rm = expand("mask/{{sample}}_repeatmasker_{n}.fa.tbl", n = N),
      rm_good = expand("mask/{{sample}}_unmaskedgood_{n}.fa", n = N),
      refgenome_unmapped = expand("refgenomefilter/{{sample}}_refgenome_unmapped_{n}.fa", n = N),
      refgenome_filtered = expand("refgenomefilter/{{sample}}_refgenome_filtered_{n}_known-host.tsv", n = N)
   output: 
      outfile = "results/{sample}_stats.json" 
   script:
      "../scripts/collect_stats.py"
