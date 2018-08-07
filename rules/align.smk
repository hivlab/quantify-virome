
## Align sequences to reference genome [11]
rule bwa_mem:
    input:
        config["ref_genome"],
        ["{sample}/10_repeatmasker_good/unmasked.{n}.fa"]
    output:
        "{sample}/11_bwa_mem/mapped.{n}.bam"
    log:
        "{sample}/logs/bwa_mem.{n}.log"
    params: "-L 100,100 -k 15"
    threads: 8
    conda:
      "envs/bwa-sam-bed.yml"
    shell:
        "(bwa mem {params} -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
