#
#bam2fastq --output {output.unmappfq} -q {output.unmapped}
## Extract unmapped reads
rule unmapped_reads:
    input:
      os.path.join(config["outdir"], dynamic("{sample}/11_bwa_mem/mapped.{n}.bam"))
    output:
      unmapped = os.path.join(config["outdir"], dynamic("{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.bam")),
      unmappfq = os.path.join(config["outdir"], dynamic("{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.fastq"))
    shell:
      """
        samtools view -b -f 4 {input} > {output.unmapped}
        bedtools bamtofastq -i {output.unmapped} -fq {output.unmappfq}
      """
