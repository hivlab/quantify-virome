
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = SAMPLES.loc[wildcards.sample, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = SAMPLES.loc[wildcards.sample, ['frac']][0]
    return frac

rule preprocess:
  input:
    sample = lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
  output:
    adapters = temp("munge/{sample}_adapters.fa"),
    merged = temp("munge/{sample}_merged.fq"),
    unmerged = temp("munge/{sample}_unmerged.fq"),
    reads = temp("munge/{sample}_reads.fq"),
    trimmed = temp("munge/{sample}_trimmed.fq")
  params:
    bbduk = "ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10 maq=10 minlen=100"
  threads: 2
  conda:
    "../envs/bbtools.yaml"
  shell:
    """
    bbmerge.sh in1={input[0]} in2={input[1]} outa={output.adapters}
    bbmerge.sh in1={input[0]} in2={input[1]} out={output.merged} outu={output.unmerged} adapters={output.adapters}
    cat {output.merged} {output.unmerged} > {output.reads}
    bbduk.sh in={output.reads} out={output.trimmed} ref={output.adapters} {params.bbduk}
    """

# Refgenome contaminant removal
rule bwa_mem_refgenome:
  input:
    reads = [rules.preprocess.output.trimmed]
  output:
    temp("mapped/{sample}_refgenome.bam")
  params:
    index = config["ref_genome"],
    extra = "-L 100,100 -k 15",
    sort = "none"
  log:
    "logs/{sample}_bwa_map_refgenome.log"
  threads: 2
  wrapper:
    "0.32.0/bio/bwa/mem"

# Concatenate merged reads and convert to fasta.
rule unmapped_fasta:
  input:
    rules.bwa_mem_refgenome.output
  output:
    fastq = temp("munge/{sample}_unmapped.fq"),
    fasta = temp("munge/{sample}_unmapped.fa")
  conda:
    "../envs/bbtools.yaml"
  shell:
    """
    reformat.sh in={input} out={output.fastq} unmappedonly primaryonly
    reformat.sh in={output.fastq} out={output.fasta} uniquenames
    """

# Calculate bam file stats.
rule refgenome_bam_stats:
    input:
      rules.bwa_mem_refgenome.output
    output:
      "stats/{sample}_refgenome_stats.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.32.0/bio/samtools/stats"

# Collect fastq stats.
rule munge_stats:
  input:
    rules.preprocess.output.trimmed, rules.unmapped_fasta.output
  output:
    "stats/{sample}_munge.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]
