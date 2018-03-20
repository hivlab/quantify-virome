
configfile: "config.yaml"

## Target rule
rule all:
    input:
        expand("fastqc/{sample}.stitched.merged_fastqc.html", sample = config["samples"]),
        expand("fastqc/{sample}.stitched.merged_fastqc.zip", sample = config["samples"]),
        expand("cdhit/{sample}.stitched.merged.prinseq.QCed.cdhit.fa", sample = config["samples"])

## Run cd-hit

rule cd_hit:
  input:
    "fasta/{sample}.stitched.merged.prinseq.fasta"
  output:
    clusters = "cdhit/{sample}.stitched.merged.prinseq.QCed.cdhit.fa",
    report = "cdhit/{sample}.stitched.merged.prinseq.QCed.cdhit.report"
  params:
    ""
  threads: 8
  shell:
    """
    cd-hit -i {input} -o {output.clusters} {params} -T {threads} > {output.report}
    """

## Convert fastq to fasta format

rule fastq2fasta:
  input: "prinseq/{sample}.stitched.merged.fq.gz"
  output: "fasta/{sample}.stitched.merged.prinseq.fasta"
  shell:
    """
    zcat | sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
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
        pair2 = "fastp/{sample}.pair2.truncated.gz",
        report = "fastp/{sample}.report.html"
    threads: 8
    shell:
        """
        fastp -i {input.R1} \
        -I {input.R2} \
        -o {output.pair1} \
        -O {output.pair2} \
        --trim-front1 5 \
        --trim-tail1 5 \
        --cut_mean_quality 25 \
        --length-required 50 \
        --low_complexity_filter \
        --complexity_threshold 8 \
        -h {output.report} \
        -w {threads}
        """
