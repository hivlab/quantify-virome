
# Run cd-hit to find and munge duplicate reads. We concatenate stitched fastq reads and convert to fasta.
rule cd_hit:
  input:
    rules.fastq_join.output
  output:
    fa = temp("munge/{sample}_merge_reads.fasta"),
    repres = temp("cdhit/{sample}_cdhit.fa"),
    report = "cdhit/{sample}_cdhit.report",
    clstr = temp("cdhit/{sample}_cdhit.fa.clstr")
  params:
    "-c 0.984 -G 0 -n 8 -d 0 -aS 0.984 -g 1 -r 1 -M 0"
  threads: 8
  conda:
    "../envs/cd-hit.yaml"
  shell:
    """
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
    cd-hit-est -i {output.fa} -o {output.repres} -T {threads} {params} > {output.report}
    """
