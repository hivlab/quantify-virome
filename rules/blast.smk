
## Blast against NT virus database
rule blastn_virus:
    input:
      db = config["virus_nt"],
      query = preprocessing("output/{sample}_refgenome_filtered_{n}_unmapped.fa")
    output:
      "output/blast/{sample}_blastn_virus_{n}.xml"
    params:
      show_gis = True,
      evalue = 1e-4,
      db_soft_mask = 100
    threads: 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blastn_virus.py"

## Filter blastn records for the cutoff value
rule parse_blastn_virus:
    input:
      rules.blastn_virus.output,
      rules.blastn_virus.input.query
    output:
      known_xml = "output/{sample}_blastn_virus_{n}_known-viral.xml",
      unmapped = "output/{sample}_blastn_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Blastx unmapped sequences against NR virus database
rule blastx_virus:
    input:
      db = config["virus_nr"],
      query = rules.parse_blastn_virus.output.unmapped
    output:
      "output/blast/{sample}_blastx_virus_{n}.xml"
    params:
      show_gis = True,
      evalue = 1e-2,
      db_soft_mask = 100
    threads: 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blastx_virus.py"

## Filter blastn records for the cutoff value
rule parse_blastx_virus:
    input:
      rules.blastx_virus.output,
      rules.blastx_virus.input.query
    output:
      known_xml = "output/{sample}_blastx_virus_{n}_known-viral.xml",
      unmapped = "output/{sample}_blastx_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-3
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Get unmasked known viral sequences
rule unmasked_viral:
  input:
    rules.parse_blastn_virus.output.known_xml,
    rules.parse_blastx_virus.output.known_xml,
    preprocessing("output/{sample}_refgenome_unmapped_{n}.fa")
  output:
    "output/{sample}_blastn_virus_{n}_known-viral.fa",
    "output/{sample}_blastx_virus_{n}_known-viral.fa"
  conda:
      "../envs/biopython.yml"
  script:
      "../scripts/unmasked_viral.py"

## Merge blast outputs
rule merge_unmasked_viral:
  input:
    rules.unmasked_viral.output
  output:
    "output/{sample}_known-viral_{n}_unmasked.fa"
  shell:
    """
    cat {input} > {output}
    """

rule bwa_mem:
    input:
        config["ref_bacteria"],
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

## Extract unmapped reads
rule unmapped_reads:
    input: rules.bwa_mem.output
    output:
      bam = temp("output/{sample}_refbacteria_unmapped_{n}.bam"),
      fq = temp("output/{sample}_refbacteria_unmapped_{n}.fq"),
      fa = "output/{sample}_refbacteria_unmapped_{n}.fa"
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
        samtools view -b -f 4 {input} > {output.bam}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

## Subset repeatmasker masked reads using unmapped ids
rule unmapped_masked:
    input: rules.unmapped_reads.output.fa, rules.repeatmasker_good.output.masked
    output:
      temp("output/{sample}_refgenome_unmapped_{n}_masked.fa")
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/unmapped_masked_ids.py"

