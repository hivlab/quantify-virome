
rule bwa_mem_refgenome:
  input:
    reads = [rules.repeatmasker_good.output.original_filt]
  output:
    temp("refgenomefilter/{sample}_refgenome_unmapped_{n}.bam")
  params:
    index = config["ref_genome"],
    extra = "-L 100,100 -k 15",
    sort = "none"
  log:
    "logs/{sample}_bwa_map_refgenome_{n}.log"
  threads: 2
  wrapper:
    "0.32.0/bio/bwa/mem"

# Extract unmapped reads and convert bam file to fastq file.
rule refgenome_unmapped_fastq:
    input:
      rules.bwa_mem_refgenome.output
    output:
      temp("refgenomefilter/{sample}_refgenome_unmapped_{n}.fq")
    params:
      "-n -f 4"
    threads: 2
    wrapper:
      "0.32.0/bio/samtools/bam2fq/interleaved"

# Convert fastq file to fasta file.
rule refgenome_unmapped:
  input:
    rules.refgenome_unmapped_fastq.output
  output:
    temp("refgenomefilter/{sample}_refgenome_unmapped_{n}.fa")
  shell:
    "cat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}"

# Calculate bam file stats
rule refgenome_bam_stats:
    input:
      rules.bwa_mem_refgenome.output
    output:
      "stats/{sample}_refgenome_stats_{n}.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.32.0/bio/samtools/stats"

## Subset repeatmasker masked reads using unmapped ids
rule refgenome_unmapped_masked:
    input: rules.refgenome_unmapped.output, rules.repeatmasker_good.output.masked_filt
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

# Collect stats
rule refgenome_unmapped_stats:
  input:
    expand(["refgenomefilter/{{sample}}_refgenome_unmapped_{n}_masked.fa", "refgenomefilter/{{sample}}_refgenome_filtered_{n}_unmapped.fa"], n = N)
  output:
    "stats/{sample}_refgenome.tsv"
  conda:
    "../envs/seqkit.yaml"
  shell:
    "seqkit stats {input} -T > {output}"
