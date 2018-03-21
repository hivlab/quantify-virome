
configfile: "config.yaml"

## Target rule
rule all:
    input:
        expand("fastqc/{sample}.stitched.merged_fastqc.html", sample = config["samples"]),
        expand("fastqc/{sample}.stitched.merged_fastqc.zip", sample = config["samples"]),
        expand("cdhit/{sample}.stitched.merged.prinseq.QCed.cdhit.fa", sample = config["samples"])

## Run FastQC -------------------------------------

rule fastqc:
    input:
        "merged/{sample}.stitched.merged.fq.gz"
    output:
        html = "fastqc/{sample}.stitched.merged_fastqc.html",
        zip = "fastqc/{sample}.stitched.merged_fastqc.zip"
    wrapper:
        "0.22.0/bio/fastqc"


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
  input: "prinseq/{sample}.stitched.merged.prinseq.fastq"
  output: "fasta/{sample}.stitched.merged.prinseq.fasta"
  shell:
    """
    sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
    rm {input}
    """

## Prinseq

rule prinseq:
    input: "merged/{sample}.stitched.merged.fq.gz"
    output:
      "prinseq/{sample}.stitched.merged.prinseq.fastq"
    params:
      settings = "-no_qual_header -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep 14 -lc_method dust -lc_threshold 8 -trim_tail_left 5 -trim_tail_right 5 -trim_ns_left 1 -trim_ns_right 1 -min_qual_mean 25 -trim_qual_left 25 -trim_qual_right 25 -min_qual_score 10",
      stub = "prinseq/{sample}.stitched.merged.prinseq"
    shell:
      """
      gzip -dc {input} | prinseq-lite.pl -fastq stdin {params.settings} -out_good {params.stub}
      rm stdin*
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
    pair1 = "adapter_removal/{sample}.pair1.truncated.gz",
    pair2 = "adapter_removal/{sample}.pair2.truncated.gz"
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

## Adapter removal --------------------------------------

rule adapter_removal:
    input:
        R1 = "data/{sample}_SE1.fastq.gz",
        R2 = "data/{sample}_SE2.fastq.gz"
    output:
        pair1 = "adapter_removal/{sample}.pair1.truncated.gz",
        pair2 = "adapter_removal/{sample}.pair2.truncated.gz",
        singletons = "adapter_removal/{sample}.singletons.truncated.gz"
    threads: 8
    shell:
        """
        AdapterRemoval \
        --threads {threads} \
        --file1 {input.R1} \
        --file2 {input.R2} \
        --output1 {output.pair1} \
        --output2 {output.pair2} \
        --singleton {output.singletons} \
        --gzip \
        --trimqualities \
        --trimns
        """
