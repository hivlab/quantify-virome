
#SAMPLES = ["SE1", "SE2"]

## Adapter removal --------------------------------------

rule adapter_removal:
    input:
        reads = "data/samples/{batch}_{sample}.fastq.gz",
        adapter_list = "data/samples/{batch}_adapter.txt"
    output:
        reads = "output/trimmed_reads/{batch}_{sample}.truncated.gz"
    shell:
        """
        AdapterRemoval \
        --file1 {input.reads} \
        --output1 {output.reads} \
        --gzip \
        --trimqualities \
        --trimns \
        --adapter-list {input.adapter_list}
        """


## Stitch paired reads --------------------------------------

rule $stitching_report_filefix:
  input:
    reads1 = "output/trimmed_reads/{batch}_SE1.truncated.gz",
    reads2 = "output/trimmed_reads/{batch}_SE2.truncated.gz"
  output:
    report = "_stitch-length-report"
  shell:
    """
    fastq-join \
    -p 5 \
    -m 10 \
    -r \${OVERLAPLENGTH}  \
    {input.reads1} \
    {input.reads2} \
    -o {batch}.%.fq
    """

