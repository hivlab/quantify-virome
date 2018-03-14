
configfile: "/gpfs/hpchome/taavi74/Projects/vs/config.yaml"

## Target rule
rule all:
    input:
        expand("output/qc/{sample}.stitched.merged_fastqc.html", sample = config["samples"]),
        expand("output/qc/{sample}.stitched.merged_fastqc.zip", sample = config["samples"])

## Run FastQC -------------------------------------

rule fastqc:
    input:
        "output/stitched/{sample}.stitched.merged.fq.gz"
    output:
        html = "output/qc/{sample}.stitched.merged_fastqc.html",
        zip = "output/qc/{sample}.stitched.merged_fastqc.zip"
    params: ""
    wrapper:
        "0.22.0/bio/fastqc"

## Merge stitched reads --------------------------------------

rule merge_reads:
  input:
    join = "output/stitched/{sample}.join.fq.gz",
    un1 = "output/stitched/{sample}.un1.fq.gz",
    un2 = "output/stitched/{sample}.un2.fq.gz"
  output:
    merged = "output/stitched/{sample}.stitched.merged.fq.gz"
  shell:
    """
    cat {input.join} {input.un1} {input.un2} > {output.merged}
    """

## Stitch paired reads --------------------------------------

rule stitch_reads:
  input:
    pair1 = "output/trimmed/{sample}.pair1.truncated.gz",
    pair2 = "output/trimmed/{sample}.pair2.truncated.gz"
  output:
    "output/stitched/{sample}.join.fq.gz",
    "output/stitched/{sample}.un1.fq.gz",
    "output/stitched/{sample}.un2.fq.gz"
  params:
    "output/stitched/{sample}.%.fq.gz"
  shell:
    """
    fastq-join -p 5 -m 10 {input.pair1} {input.pair2} -o {params}
    """

## Adapter removal --------------------------------------

rule adapter_removal:
    input:
        R1 = "raw/{sample}_SE1.fastq.gz",
        R2 = "raw/{sample}_SE2.fastq.gz"
    output:
        pair1 = "output/trimmed/{sample}.pair1.truncated.gz",
        pair2 = "output/trimmed/{sample}.pair2.truncated.gz",
        singletons = "output/trimmed/{sample}.singletons.truncated.gz"
    shell:
        """
        AdapterRemoval \
        --file1 {input.R1} \
        --file2 {input.R2} \
        --output1 {output.pair1} \
        --output2 {output.pair2} \
        --singleton {output.singletons} \
        --gzip \
        --trimqualities \
        --trimns
        """
