
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = SAMPLES.loc[wildcards.sample, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = SAMPLES.loc[wildcards.sample, ['frac']][0]
    return frac

# Adapter trimming and quality filtering.
rule fastp:
  input:
    sample = lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
  output:
    trimmed = [temp("munge/{sample}_trimmed1.fq"), temp("munge/{sample}_trimmed2.fq")],
    json = "stats/{sample}_fastp.json",
    html = "stats/{sample}_fastp.html"
  params:
    extra = "--trim_front1 5 --trim_tail1 5 --length_required 50 --low_complexity_filter --complexity_threshold 8"
  threads: 2
  log:
    "logs/{sample}_fastp.log"
  wrapper:
    "0.34.0/bio/fastp"

# Refgenome contaminant removal
rule bwa_mem_refgenome:
  input:
    reads = [rules.fastp.output.trimmed]
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

# Extract unmapped reads and convert bam file to fastq file.
rule refgenome_unmapped_fastq:
  input:
    rules.bwa_mem_refgenome.output
  output:
    temp("munge/{sample}_unmapped1.fq"),
    temp("munge/{sample}_unmapped2.fq")
  params:
    sort = "",
    bam2fq = "-f 4"
  threads: 2
  wrapper:
    "0.32.0/bio/samtools/bam2fq/separate"

# Subsample much bigger runs
# based on precalculated fractions in samples.tsv.
rule sample:
  input:
    rules.refgenome_unmapped_fastq.output
  output:
    temp("munge/{sample}_subsample1.fq"),
    temp("munge/{sample}_subsample2.fq")
  params:
    frac = get_frac,
    seed = config["seed"]
  wrapper:
    config["wrappers"]["sample"]

# Stitch paired reads.
rule fastq_join:
  input:
    rules.sample.output
  output:
    temp("munge/{sample}_un1.fq"),
    temp("munge/{sample}_un2.fq"),
    temp("munge/{sample}_join.fq")
  params:
    options = "-p 5 -m 10 -r stats/{sample}_stitchlength.report"
  log:
    "logs/{sample}_fastq_join.log"
  wrapper:
    config["wrappers"]["fastq_join"]

# Calculate bam file stats.
rule refgenome_bam_stats:
    input:
      rules.bwa_mem_refgenome.output
    output:
      "stats/{sample}_refgenome_mapping.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.32.0/bio/samtools/stats"

# Collect fastq stats.
rule munge_stats:
  input:
    rules.fastp.output.trimmed, rules.refgenome_unmapped_fastq.output, rules.sample.output, rules.fastq_join.output
  output:
    "stats/{sample}_munge.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]
