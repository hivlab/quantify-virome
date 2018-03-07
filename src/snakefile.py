
## Adapter removal --------------------------------------

rule adapter_removal:
    input:
        reads = "data/samples/{sample}_{read}.fastq.gz",
        adapter_list = "data/samples/{sample}_adapter.txt"
    output:
        reads = "output/trimmed_reads/{sample}_{read}.truncated.gz"
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
    reads1 = "output/trimmed_reads/{sample}_SE1.truncated.gz",
    reads2 = "output/trimmed_reads/{sample}_SE2.truncated.gz"
  output:
    report = "output/stitched_reads/{sample}.stitch-length-report",
    dir = "output/stitched_reads/"
  shell:
    """
    fastq-join \
    -p 5 \
    -m 10 \
    -r {output.report}  \
    {input.reads1} \
    {input.reads2} \
    -o {output.dir}/{wildcards.batch}.%.fq.gz
    """

