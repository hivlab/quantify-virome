
configfile: "config.yaml"

## Target rule
rule all:
    input:
        expand("cdhit/{sample}.stitched.merged.cdhit.fa", sample = config["samples"]),
        expand("cdhit/{sample}.stitched.merged.cdhit.report", sample = config["samples"])

## Run cd-hit

rule cd_hit:
  input:
    "fasta/{sample}.stitched.merged.fasta"
  output:
    clusters = "cdhit/{sample}.stitched.merged.cdhit.fa",
    report = "cdhit/{sample}.stitched.merged.cdhit.report"
  params:
    "-M 10000"
  threads: 8
  shell:
    """
    cd-hit -i {input} -o {output.clusters} {params} -T {threads} > {output.report}
    """

## Convert fastq to fasta format

rule fastq2fasta:
  input: "merged/{sample}.stitched.merged.fq.gz"
  output: "fasta/{sample}.stitched.merged.fasta"
  shell:
    """
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
    rm {input}
    """


## Merge stitched reads --------------------------------------

rule merge_reads:
  input:
    join = "stitched/{sample}.join.fq.gz",
    un1 = "stitched/{sample}.un1.fq.gz",
    un2 = "stitched/{sample}.un2.fq.gz"
  output:
    merged = "merged/{sample}.stitched.merged.fq.gz"
  shell:
    """
    cat {input.join} {input.un1} {input.un2} > {output.merged}
    """

## Stitch paired reads --------------------------------------

rule stitch_reads:
  input:
    pair1 = "fastp/{sample}.pair1.truncated.gz",
    pair2 = "fastp/{sample}.pair2.truncated.gz"
  output:
    "stitched/{sample}.join.fq.gz",
    "stitched/{sample}.un1.fq.gz",
    "stitched/{sample}.un2.fq.gz"
  params:
    maximum_difference = 5,
    minimum_overlap = 10,
    template = "stitched/{sample}.%.fq.gz"
  shell:
    """
    fastq-join \
    -p {params.maximum_difference} \
    -m {params.minimum_overlap} \
    {input.pair1} \
    {input.pair2} \
    -o {params.template}
    """

## All-in-one preprocessing for FastQ files -----------------------------------
# Adapter trimming is enabled by default
# Quality filtering is enabled by default

rule fastp:
    input:
        R1 = "data/{sample}_SE1.fastq.gz",
        R2 = "data/{sample}_SE2.fastq.gz"
    output:
        pair1 = "fastp/{sample}.pair1.truncated.gz",
        pair2 = "fastp/{sample}.pair2.truncated.gz"
    params:
        report = "fastp/{sample}.report.html"
    threads: 8
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.pair1} -O {output.pair2} -f 5 -t 5 -M 25 -l 50 -y -Y 8 -h {params.report} -w {threads}
        """
