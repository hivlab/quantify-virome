
import os.path

configfile: "config.yaml"

## Target rule
rule all:
    input:
        expand(os.path.join(config["outdir"], "cdhit/{sample}.stitched.merged.cdhit.fa"), sample = config["samples"]),
        expand(os.path.join(config["outdir"], "cdhit/{sample}.stitched.merged.cdhit.report"), sample = config["samples"])

## Run cd-hit

rule cd_hit:
  input:
    os.path.join(config["outdir"], "fasta/{sample}.stitched.merged.fasta")
  output:
    clusters = os.path.join(config["outdir"], "cdhit/{sample}.stitched.merged.cdhit.fa"),
    report = os.path.join(config["outdir"], "cdhit/{sample}.stitched.merged.cdhit.report")
  resources:
    mem = 8
  params:
    "-M 8000"
  threads:
    20
  shell:
    """
    cd-hit -i {input} -o {output.clusters} {params} -T {threads} {params} > {output.report}
    """

## Convert fastq to fasta format

rule fastq2fasta:
  input: os.path.join(config["outdir"], "merged/{sample}.stitched.merged.fq.gz")
  output: os.path.join(config["outdir"], "fasta/{sample}.stitched.merged.fasta")
  shell:
    """
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
    rm {input}
    """


## Merge stitched reads --------------------------------------

rule merge_reads:
  input:
    join = os.path.join(config["outdir"], "stitched/{sample}.join.fq.gz"),
    un1 = os.path.join(config["outdir"], "stitched/{sample}.un1.fq.gz"),
    un2 = os.path.join(config["outdir"], "stitched/{sample}.un2.fq.gz")
  output:
    merged = os.path.join(config["outdir"], "merged/{sample}.stitched.merged.fq.gz")
  shell:
    """
    cat {input.join} {input.un1} {input.un2} > {output.merged}
    """

## Stitch paired reads --------------------------------------

rule fastq_join:
  input:
    pair1 = os.path.join(config["outdir"], "fastp/{sample}.pair1.truncated.gz"),
    pair2 = os.path.join(config["outdir"], "fastp/{sample}.pair2.truncated.gz")
  output:
    os.path.join(config["outdir"], "stitched/{sample}.join.fq.gz"),
    os.path.join(config["outdir"], "stitched/{sample}.un1.fq.gz"),
    os.path.join(config["outdir"], "stitched/{sample}.un2.fq.gz")
  params:
    maximum_difference = 5,
    minimum_overlap = 10,
    template = os.path.join(config["outdir"], "stitched/{sample}.%.fq.gz")
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
# Replaces AdapteRemoval, prinseq and fastqc

rule fastp:
    input:
        R1 = os.path.join(config["datadir"], "{sample}_SE1.fastq.gz"),
        R2 = os.path.join(config["datadir"], "{sample}_SE2.fastq.gz")
    output:
        pair1 = os.path.join(config["outdir"], "fastp/{sample}.pair1.truncated.gz"),
        pair2 = os.path.join(config["outdir"], "fastp/{sample}.pair2.truncated.gz")
    params:
        report = "fastp/{sample}.report.html"
    threads: 8
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.pair1} -O {output.pair2} -f 5 -t 5 -M 25 -l 50 -y -Y 8 -h {params.report} -w {threads}
        """
