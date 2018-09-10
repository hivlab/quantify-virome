
def get_fastq(wildcards, read_pair = 'fq1'):
 return samples.loc[wildcards.sample, [read_pair]].dropna()[0]

## Preprocessing for fastq files
# Adapter trimming
# Quality filtering
# Replaces AdapteRemoval, prinseq and fastqc
rule fastp:
    input:
      fq1 = lambda wildcards: get_fastq(wildcards, 'fq1'),
      fq2 = lambda wildcards: get_fastq(wildcards, 'fq2')
    output:
      pair1 = "output/{sample}_pair1_trimmed.gz",
      pair2 = "output/{sample}_pair2_trimmed.gz",
      html = "output/logs/{sample}_fastp_report.html",
      json = "output/logs/{sample}_fastp_report.json"
    params:
      "-f 5 -t 5 -l 50 -y -Y 8"
    threads: 8
    conda:
      "../envs/fastp.yml"
    log: "output/logs/{sample}_fastp.log"
    shell:
      """
      fastp -i {input.fq1} -I {input.fq2} -o {output.pair1} -O {output.pair2} {params} -h {output.html} -j {output.json} -w {threads} > {log} 2>&1
      """

## Stitch paired reads
rule fastq_join:
  input: rules.fastp.output
  output:
    "output/{sample}_join.fq.gz",
    "output/{sample}_un1.fq.gz",
    "output/{sample}_un2.fq.gz"
  params:
    config["fastq-join"]["maximum_difference"],
    config["fastq-join"]["minimum_overlap"],
    "output/{sample}_%.fq.gz"
  conda:
    "../envs/fastq-join.yml"
  log: "output/logs/{sample}_fastq_join.log"
  shell:
    """
    fastq-join \
    -p {params[0]} \
    -m {params[1]} \
    {input[0]} \
    {input[1]} \
    -o {params[2]} > {log} 2>&1
    """

## Merge stitched reads
rule merge_reads:
  input: rules.fastq_join.output
  output:
    fq = "output/{sample}_merge_reads.fq.gz",
    fa = "output/{sample}_merge_reads.fasta"
  shell:
    """
    cat {input} > {output.fq}
    zcat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
    """
