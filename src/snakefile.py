
configfile: "config.yaml"

## Adapter removal --------------------------------------

rule adapter_removal:
    input:
        reads1 = "data/samples/{sample}_SE1.fastq.gz",
        reads2 = "data/samples/{sample}_SE2.fastq.gz"#,
        # adapter_list = "data/samples/{sample}_adapter.txt"
    output:
        pair1 = "output/trimmed_reads/{sample}.pair1.truncated.gz",
        pair2 = "output/trimmed_reads/{sample}.pair2.truncated.gz",
        singleton = "output/trimmed_reads/{sample}.singletons.truncated.gz"
    shell:
        """
        AdapterRemoval \
        --file1 {input.reads1} \
        --file2 {input.reads2} \
        --output1 {output.pair1} \
        --output2 {output.pair2} \
        --singleton {output.singleton} \
        --gzip \
        --trimqualities \
        --trimns \
        --trim5p 9 \
        --trim3p 9
        """

## Stitch paired reads --------------------------------------

rule stitch_reads:
  input:
    pair1 = "output/trimmed_reads/{sample}.pair1.truncated.gz",
    pair2 = "output/trimmed_reads/{sample}.pair2.truncated.gz"
  output:
    report = "output/stitched_reads/{sample}.stitch-length-report",
    out = "output/stitched_reads/{sample}.%.fq.gz"
  shell:
    """
    fastq-join \
    -p 5 \
    -m 10 \
    -r {output.report} \
    {input.pair1} \
    {input.pair2} \
    -o {output.out}
    """

## Merge stitched reads --------------------------------------

rule merge_reads:
  input:
    join = "output/stitched_reads/{sample}.join.fq.gz",
    un1 = "output/stitched_reads/{sample}.un1.fq.gz",
    un2 = "output/stitched_reads/{sample}.un2.fq.gz"
  output:
    merged = "output/stitched_reads/{sample}.stitched.merged.fq.gz"
  shell:
    """
    cat {input.join} {input.un1} {input.un2} > {output.merged}
    """

## Run FastQC -------------------------------------

rule fastqc:
    input:
        "output/stitched_reads/{sample}.stitched.merged.fq.gz"
    output:
        "output/qc/{sample}.stitched.merged_fastqc.html",
        "output/qc/{sample}.stitched.merged_fastqc.zip"
    params: ""
    shell:
        "fastqc --outdir output/qc/ {input}"
