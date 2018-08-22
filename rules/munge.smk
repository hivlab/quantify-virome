
def get_fastq(wildcards, read_pair = 'fq1'):
 return samples.loc[wildcards.sample, [read_pair]].dropna()[0]

## All-in-one preprocessing for FastQ files [1,3]
# Adapter trimming is enabled by default
# Quality filtering is enabled by default
# Replaces AdapteRemoval, prinseq and fastqc
rule fastp:
    input:
      fq1 = lambda wildcards: get_fastq(wildcards, 'fq1'),
      fq2 = lambda wildcards: get_fastq(wildcards, 'fq2')
    output:
      pair1 = temp("output/01_fastp/{sample}_pair1_truncated.gz"),
      pair2 = temp("output/01_fastp/{sample}_pair2_truncated.gz"),
      html = "output/01_fastp/{sample}_report.html",
      json = "output/01_fastp/{sample}_report.json"
    params:
      options = "-f 5 -t 5 -l 50 -y -Y 8"
    threads: 8
    conda:
      "../envs/fastp.yml"
    log: "output/logs/{sample}_fastp.log"
    shell:
      """
      fastp -i {input.fq1} -I {input.fq2} -o {output.pair1} -O {output.pair2} {params.options} -h {output.html} -j {output.json} -w {threads} > {log} 2>&1
      """

## Stitch paired reads [2]
rule fastq_join:
  input: rules.fastp.output
  output:
    temp("output/02_stitched/{sample}_join.fq.gz"),
    temp("output/02_stitched/{sample}_un1.fq.gz"),
    temp("output/02_stitched/{sample}_un2.fq.gz")
  params:
    config["fastq-join"]["maximum_difference"],
    config["fastq-join"]["minimum_overlap"],
    template = "output/02_stitched/{sample}_%.fq.gz"
  conda:
    "../envs/fastq-join.yml"
  log: "{sample}/logs/02_fastq_join.log"
  shell:
    """
    fastq-join \
    -p {params[0]} \
    -m {params[1]} \
    {input[0]} \
    {input[1]} \
    -o {params.template} > {log} 2>&1
    """

## Merge stitched reads [3]
rule merge_reads:
  input: rules.fastq_join.output
  output:
    temp("output/03_merged/{sample}_stitched_merged.fq.gz")
  shell:
    """
    cat {input} > {output}
    """

## Convert fastq to fasta format [4]
rule fastq2fasta:
  input: rules.merge_reads.output
  output:
    temp("sample/04_fasta/{sample}_stitched_merged.fasta")
  shell:
    """
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
    """

## Run cd-hit to find and munge duplicate reads [5]
rule cd_hit:
  input: rules.fastq2fasta.output
  output:
    clusters = "output/05_cdhit/{sample}_merged_cdhit.fa",
    report = "output/05_cdhit/{sample}_merged_cdhit.report"
  threads: 8
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 2000
  conda:
    "../envs/cd-hit.yml"
  shell:
    """
    cd-hit-est -i {input} -o {output.clusters} -T {threads} -c 0.984 -G 0 -n 8 -d 0 -aS 0.984 -g 1 -r 1 -M 0 > {output.report}
    """
