
#SAMPLES = ["SE1", "SE2"]

##--------------------------------------##
## 1. Adapter removal                   ##
##--------------------------------------##

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
