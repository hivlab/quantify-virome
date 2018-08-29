
## Align sequences to reference genome [11]
rule bwa_mem:
    input:
        config["ref_genome"],
        ["output/{sample}_unmaskedgood_{n}.fa"]
    output:
        "output/bwa_mem/{sample}_mapped_{n}.bam"
    log:
        "output/logs/{sample}_bwa_mem_{n}.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
        "(bwa mem -L 100,100 -k 15 -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

## Extract unmapped reads [12a]
rule unmapped_reads:
    input: rules.bwa_mem.output
    output:
      bam = temp("output/{sample}_refgenome_unmapped_{n}.bam"),
      fq = temp("output/{sample}_refgenome_unmapped_{n}.fq"),
      fa = temp("output/{sample}_refgenome_unmapped_{n}.fa")
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
        samtools view -b -f 4 {input} > {output.bam}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

## Subset repeatmasker masked reads using unmapped ids [12b]
rule unmapped_masked:
    input: rules.unmapped_reads.output.fa, rules.repeatmasker_good.output.masked
    output:
      temp("output/{sample}_refgenome_unmapped_{n}_masked.fa")
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/unmapped_masked_ids.py"
