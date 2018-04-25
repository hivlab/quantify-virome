
## Extract unmapped reads [12]
rule unmapped_reads:
    input:
      os.path.join(config["outdir"], dynamic("{sample}/11_bwa_mem/mapped.{n}.bam"))
    output:
      bam = os.path.join(config["outdir"], dynamic("{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.bam")),
      fq = os.path.join(config["outdir"], dynamic("{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.fq")),
      fa = os.path.join(config["outdir"], dynamic("{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.fa"))
    shell:
      """
        samtools view -b -f 4 {input} > {output.bam}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

rule subset_unmapped:
    input:
      os.path.join(config["outdir"], dynamic("{sample}/12_unmapped_reads/RefGenome_unmapped.{n}.fa")),
      os.path.join(config["outdir"], dynamic("{sample}/10_repeatmasker_good/masked.{n}.fa"))
    output:
      os.path.join(config["outdir"], dynamic("{sample}/13_unmapped_masked/RefGenome_unmapped.{n}.masked.fa"))
    script:
      "../scripts/unmapped_masked_ids.py"
