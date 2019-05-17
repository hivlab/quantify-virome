
# Run cd-hit to find and munge duplicate reads.
rule cd_hit:
  input:
    rules.unmapped_fasta.output.fasta
  output:
    repres = temp("cdhit/{sample}_cdhit.fa"),
    clstr = temp("cdhit/{sample}_cdhit.fa.clstr")
  params:
    "-c 0.984 -G 0 -n 10 -d 0 -aS 0.984 -r 1 -M 0"
  threads: 2
  log:
    "logs/{sample}_cdhit.log"
  conda:
    "../envs/cd-hit.yaml"
  shell:
    """
    cd-hit-est -i {input} -o {output.repres} -T {threads} {params} > {log}
    """
