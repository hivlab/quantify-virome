
## Align sequences to reference genome and extract unmapped reads
rule refgenome_unmapped:
    input:
        config["ref_genome"],
        [rules.repeatmasker_good.output.original_filt]
    output:
      bam = temp("refgenomefilter/{sample}_refgenome_unmapped_{n}.bam"),
      fq = temp("refgenomefilter/{sample}_refgenome_unmapped_{n}.fq"),
      fa = "refgenomefilter/{sample}_refgenome_unmapped_{n}.fa"
    log:
        "logs/{sample}_bwa_map_refgenome_{n}.log"
    threads: 2
    conda:
      "../envs/bwa-sam-bed.yaml"
    shell:
      """
        (bwa mem -L 100,100 -k 15 -t {threads} {input} | samtools view -b -S -f 4 - > {output.bam}) 2> {log}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

## Subset repeatmasker masked reads using unmapped ids
rule refgenome_unmapped_masked:
    input: rules.refgenome_unmapped.output.fa, rules.repeatmasker_good.output.masked_filt
    output:
      temp("refgenomefilter/{sample}_refgenome_unmapped_{n}_masked.fa")
    conda:
      "../envs/biopython.yaml"
    script:
      "../scripts/unmapped_masked_ids.py"

## MegaBlast against reference genome to remove host sequences
rule megablast_refgenome:
    input:
      query = rules.refgenome_unmapped_masked.output
    output:
      out = "refgenomefilter/{sample}_megablast_{n}.tsv"
    params:
      db = config["ref_genome"],
      task = "megablast",
      perc_identity = config["megablast_ref_genome"]["perc_identity"],
      evalue = config["megablast_ref_genome"]["evalue"],
      word_size = config["megablast_ref_genome"]["word_size"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = "'6 qseqid sgi pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    wrapper:
      config["wrappers"]["blast"]

## Filter megablast records for the cutoff value
rule parse_megablast:
    input:
      blast_result = rules.megablast_refgenome.output.out,
      query = rules.refgenome_unmapped_masked.output
    output:
      mapped = "refgenomefilter/{sample}_refgenome_filtered_{n}_known-host.tsv",
      unmapped = temp("refgenomefilter/{sample}_refgenome_filtered_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]
