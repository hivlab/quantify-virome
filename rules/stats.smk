
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

if config["zenodo"]["deposition_id"]:
   localrules: all, upload_stats
   rule upload_stats:
      input: rules.collect_stats.output
      output: ZEN.remote(expand("{deposition_id}/files/{{sample}}_stats.json", deposition_id = config["zenodo"]["deposition_id"]))
      shell: "cp {input} {output}"
