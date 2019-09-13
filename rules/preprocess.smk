
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = RUNS.loc[wildcards.run, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = RUNS.loc[wildcards.run, ['frac']][0]
    return frac

rule preprocess:
  input:
    lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
  output:
    adapters = temp("preprocess/{run}_adapters.fa"),
    merged = temp("preprocess/{run}_merged.fq"),
    unmerged = temp("preprocess/{run}_unmerged.fq"),
    reads = temp("preprocess/{run}_reads.fq"),
    trimmed = temp("preprocess/{run}_trimmed.fq"),
    sampled = temp("preprocess/{run}_sample.fq")
  params:
    bbduk = "qtrim=r trimq=10 maq=10 minlen=100",
    frac = lambda wildcards: get_frac(wildcards),
    seed = config["seed"]
  threads: 2
  wrapper:
    "https://raw.githubusercontent.com/avilab/virome-wrappers/master/preprocess"

# Map reads to Refgenome.
rule bwa_mem_refgenome:
  input:
    reads = [rules.preprocess.output.sampled]
  output:
    temp("mapped/{run}_refgenome.bam")
  params:
    index = REF_GENOME,
    extra = "-L 100,100 -k 15",
    sort = "none"
  log:
    "logs/{run}_bwa_map_refgenome.log"
  threads: 2
  wrapper:
    "0.32.0/bio/bwa/mem"

# Extract unmapped reads and convert to fasta.
rule unmapped_refgenome:
  input:
    rules.bwa_mem_refgenome.output
  output:
    fastq = temp("preprocess/{run}_unmapped.fq"),
    fasta = temp("preprocess/{run}_unmapped.fa")
  params:
    reformat_fasta_extra = "uniquenames",
    extra = "-Xmx48000m"
  wrapper:
    BWA_UNMAPPED

# Run cd-hit to find and cluster duplicate reads.
rule cd_hit:
  input:
    rules.unmapped_refgenome.output.fasta
  output:
    repres = temp("cdhit/{run}_cdhit.fa"),
    clstr = temp("cdhit/{run}_cdhit.fa.clstr")
  params:
    extra = "-c 0.984 -G 0 -n 10 -d 0 -aS 0.984 -r 1 -M 0"
  threads: 2
  log:
    "logs/{run}_cdhit.log"
  wrapper:
    "https://raw.githubusercontent.com/avilab/virome-wrappers/master/cdhit"

# Tantan mask of low complexity DNA sequences
rule tantan:
  input:
    rules.cd_hit.output.repres
  output:
    temp("mask/{run}_tantan.fasta")
  params:
    extra = "-x N" # mask low complexity using N
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/7e681180a5607f20594b3070f8eced7ccd245a89/bio/tantan"

# Filter tantan output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
  input:
    masked = rules.tantan.output
  output:
    masked_filt = temp("mask/{run}_tantangood.fasta")
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    LN_FILTER

# Split reads to smaller chunks for Repeatmasker
rule split_fasta:
  input:
    rules.tantan_good.output
  output:
    temp(expand("mask/{{run}}_repeatmasker_{n}.fa", n = N))
  params:
    config["split_fasta"]["n_files"]
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/7e681180a5607f20594b3070f8eced7ccd245a89/bio/split-fasta"

# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
  input:
    fa = "mask/{run}_repeatmasker_{n}.fa"
  output:
    masked = temp("mask/{run}_repeatmasker_{n}.fa.masked"),
    out = temp("mask/{run}_repeatmasker_{n}.fa.out"),
    cat = temp("mask/{run}_repeatmasker_{n}.fa.cat"),
    tbl = "mask/{run}_repeatmasker_{n}.fa.tbl"
  params:
    outdir = "mask"
  threads: 2
  singularity:
    "shub://tpall/repeatmasker-singularity"
  shell:
    """
    RepeatMasker -qq -pa {threads} {input.fa} -dir {params.outdir}
    # Keep consistent outputs
    # - Generate missing .tbl file when no repetitive seqs were detected 
    if head -n 1 {output.out} | grep -q "There were no repetitive sequences detected"
    then 
      ln -sr {input.fa} {output.masked} && touch {output.tbl}
    fi
    # - Unzip cat.gz file that is created if totseqlen > 10000000
    if [[ -f mask/{wildcards.run}_repeatmasker_{wildcards.n}.fa.cat.gz ]]
    then
      gzip -d mask/{wildcards.run}_repeatmasker_{wildcards.n}.fa.cat.gz
    fi
    """

# Filter repeatmasker output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
  input:
    masked = rules.repeatmasker.output.masked,
    original = rules.repeatmasker.input.fa
  output:
    masked_filt = temp("mask/{run}_repmaskedgood_{n}.fa"),
    original_filt = temp("mask/{run}_unmaskedgood_{n}.fa")
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    LN_FILTER

# MegaBlast against reference genome to remove host sequences
rule megablast_refgenome:
    input:
      query = rules.repeatmasker_good.output.masked_filt
    output:
      out = temp("blast/{run}_megablast_{n}.tsv")
    params:
      program = "blastn",
      db = REF_GENOME,
      task = "megablast",
      perc_identity = 85,
      evalue = 1e-10,
      word_size = 16,
      max_hsps = 1,
      outfmt = "'6 qseqid sseqid pident length evalue'"
    threads: 2
    wrapper:
      BLAST_QUERY

# Filter megablast records for the cutoff value
rule parse_megablast_refgenome:
    input:
      blast_result = rules.megablast_refgenome.output.out,
      query = rules.repeatmasker_good.output.masked_filt
    output:
      mapped = temp("blast/{run}_refgenome-megablast_{n}_known-host.tsv"),
      unmapped = temp("blast/{run}_refgenome-megablast_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      PARSE_BLAST

# Collect stats from preprocess outputs.
rule preprocess_stats:
  input:
    rules.preprocess.output.trimmed,
    rules.unmapped_refgenome.output,
    expand("blast/{{run}}_refgenome-megablast_{n}_unmapped.fa", n = N),
    rules.cd_hit.output.repres,
    rules.tantan.output,
    rules.tantan_good.output,
    expand(["mask/{{run}}_repmaskedgood_{n}.fa", "mask/{{run}}_unmaskedgood_{n}.fa"], n = N)
  output:
    "stats/{run}_preprocess.tsv"
  params:
    extra = "-T"
  wrapper:
    SEQ_STATS

# Refgenome mapping stats.
rule refgenome_bam_stats:
    input:
      rules.bwa_mem_refgenome.output
    output:
      "stats/{run}_refgenome-stats.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.32.0/bio/samtools/stats"
