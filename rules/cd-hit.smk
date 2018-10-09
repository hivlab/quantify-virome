
## Run cd-hit to find and munge duplicate reads
rule cd_hit:
  input: rules.merge_reads.output.fa
  output:
    clusters = "cdhit/{sample}_cdhit.fa",
    report = "cdhit/{sample}_cdhit.report",
    clstr = "cdhit/{sample}_cdhit.fa.clstr"
  params:
    "-c 0.984 -G 0 -n 8 -d 0 -aS 0.984 -g 1 -r 1 -M 0"
  threads: 8
  conda:
    "../envs/cd-hit.yml"
  shell:
    """
    cd-hit-est -i {input} -o {output.clusters} -T {threads} {params} > {output.report}
    """
