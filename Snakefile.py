
configfile: "config.yaml"

SAMPLES = config["samples"]
RAW = config["rawseqs"]
OUT = config["output"]

## Target rule
rule all:
    input:
        expand("{out}/qc/{sample}.stitched.merged_fastqc.html", out = OUT, sample = SAMPLES),
        expand("{out}/qc/{sample}.stitched.merged_fastqc.zip", out = OUT, sample = SAMPLES)

## Adapter removal --------------------------------------

rule adapter_removal:
    input:
        R1 = expand("{raw}/{sample}_SE1.fastq.gz", raw = RAW, sample = SAMPLES),
        R2 = expand("{raw}/{sample}_SE2.fastq.gz", raw = RAW, sample = SAMPLES)
    output:
        pair1 = expand("{out}/trimmed/{sample}.pair1.truncated.gz", out = OUT, sample = SAMPLES),
        pair2 = expand("{out}/trimmed/{sample}.pair2.truncated.gz", out = OUT, sample = SAMPLES),
        singletons = expand("{out}/trimmed/{sample}.singletons.truncated.gz", out = OUT, sample = SAMPLES)
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

## Stitch paired reads --------------------------------------

rule stitch_reads:
  input:
    pair1 = expand("{out}/trimmed/{sample}.pair1.truncated.gz", out = OUT, sample = SAMPLES),
    pair2 = expand("{out}/trimmed/{sample}.pair2.truncated.gz", out = OUT, sample = SAMPLES)
  output:
    template = expand("{sample}.%.fq.gz", out = OUT, sample = SAMPLES)
  shell:
    """
    fastq-join \
    -p 5 \
    -m 10 \
    {input.pair1} \
    {input.pair2} \
    -o {output.template}
    """

## Merge stitched reads --------------------------------------

rule merge_reads:
  input:
    join = expand("{out}/stitched/{sample}.join.fq.gz", out = OUT, sample = SAMPLES),
    un1 = expand("{out}/stitched/{sample}.un1.fq.gz", out = OUT, sample = SAMPLES),
    un2 = expand("{out}/stitched/{sample}.un2.fq.gz", out = OUT, sample = SAMPLES)
  output:
    merged = expand("{out}/stitched/{sample}.stitched.merged.fq.gz", out = OUT, sample = SAMPLES)
  shell:
    """
    cat {input.join} {input.un1} {input.un2} > {output.merged}
    """

## Run FastQC -------------------------------------

rule fastqc:
    input:
        merged = expand("{out}/stitched/{sample}.stitched.merged.fq.gz", out = OUT, sample = SAMPLES)
    output:
        expand("{out}/qc/{sample}.stitched.merged_fastqc.html", out = OUT, sample = SAMPLES),
        expand("{out}/qc/{sample}.stitched.merged_fastqc.zip", out = OUT, sample = SAMPLES)
    params: ""
    shell:
        "fastqc {input}"
