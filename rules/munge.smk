
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = SAMPLES.loc[wildcards.sample, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = SAMPLES.loc[wildcards.sample, ['frac']][0]
    return frac

# Imports local or remote fastq(.gz) files. Downsamples runs based on user-provided fractions in samples.tsv file.
rule sample:
  input:
    lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq
  output:
    temp("munge/{sample}_read1.fq.gz"),
    temp("munge/{sample}_read2.fq.gz")
  params:
    frac = get_frac,
    seed = config["seed"]
  wrapper:
    config["wrappers"]["sample"]

# Adapter trimming and quality filtering.
rule fastp:
  input:
    rules.sample.output
  output:
    temp("munge/{sample}_read1_trimmed.fq.gz"),
    temp("munge/{sample}_read2_trimmed.fq.gz")
  params:
    options = "--trim_front1 5 --trim_tail1 5 --length_required 50 --low_complexity_filter --complexity_threshold 8",
    json = "stats/{sample}_fastp.json"
  threads: 2
  log:
    "logs/{sample}_fastp.log"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/c27e7d713115cec135cf0123601f4b9e7799ce20/bio/fastp/wrapper.py"

# Stitch paired reads.
rule fastq_join:
  input:
    rules.fastp.output
  output:
    temp("munge/{sample}_un1.fq.gz"),
    temp("munge/{sample}_un2.fq.gz"),
    temp("munge/{sample}_join.fq.gz")
  params:
    options = "-p 5 -m 10"
  log:
    "logs/{sample}_fastq_join.log"
  wrapper:
    config["wrappers"]["fastq_join"]

# Collect fastq stats
rule munge_stats:
  input:
    rules.sample.output, rules.fastp.output, rules.fastq_join.output
  output:
    "stats/{sample}_munge.tsv"
  params:
    extra = "-T"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/dfff20d4f55ed7b9e52afa34f57a4556e295680f/bio/seqkit/stats"
