
## Adapter removal --------------------------------------

rule adapter_removal:
    input:
        reads1 = "data/samples/{sample}_SE1.fastq.gz",
        reads2 = "data/samples/{sample}_SE2.fastq.gz",
        adapter_list = "data/samples/{sample}_adapter.txt"
    output:
        pair1 = "output/trimmed_reads/{sample}.pair1.truncated.gz",
        pair2 = "output/trimmed_reads/{sample}.pair2.truncated.gz",
        singleton = "output/trimmed_reads/{sample}.singletons.truncated.gz",
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
        --adapter-list {input.adapter_list}
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
