
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
    adapters = temp("preprocess/{sample}_adapters.fa"),
    merged = temp("preprocess/{sample}_merged.fq"),
    unmerged = temp("preprocess/{sample}_unmerged.fq"),
    reads = temp("preprocess/{sample}_reads.fq"),
    trimmed = temp("preprocess/{sample}_trimmed.fq")
  params:
    bbduk = "ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10 maq=10 minlen=100"
  threads: 2
  singularity:
    "docker://bryce911/bbtools"
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
    fastq = temp("preprocess/{sample}_unmapped.fq"),
    fasta = temp("preprocess/{sample}_unmapped.fa")
  singularity:
    "docker://bryce911/bbtools"
  shell:
    """
    reformat.sh in={input} out={output.fastq} unmappedonly primaryonly
    reformat.sh in={output.fastq} out={output.fasta} uniquenames
    """

# Run cd-hit to find and cluster duplicate reads.
rule cd_hit:
  input:
    rules.unmapped_fasta.output.fasta
  output:
    repres = temp("cdhit/{sample}_cdhit.fa"),
    clstr = temp("cdhit/{sample}_cdhit.fa.clstr")
  params:
    "-c 0.984 -G 0 -n 10 -d 0 -aS 0.984 -r 1 -M 0"
  threads: 2
  log:
    "logs/{sample}_cdhit.log"
  singularity:
    "docker://quay.io/biocontainers/cd-hit:4.8.1--hdbcaa40_2"
  shell:
    """
    cd-hit-est -i {input} -o {output.repres} -T {threads} {params} > {log}
    """

# Tantan mask of low complexity DNA sequences
rule tantan:
  input:
    rules.cd_hit.output.repres
  output:
    temp("mask/{sample}_tantan.fasta")
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
    masked_filt = temp("mask/{sample}_tantangood.fasta")
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/filter/masked"

# Split reads to smaller chunks for Repeatmasker
rule split_fasta:
  input:
    rules.tantan_good.output
  output:
    temp(expand("mask/{{sample}}_repeatmasker_{n}.fa", n = N))
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
    fa = "mask/{sample}_repeatmasker_{n}.fa"
  output:
    masked = temp("mask/{sample}_repeatmasker_{n}.fa.masked"),
    out = temp("mask/{sample}_repeatmasker_{n}.fa.out"),
    tbl = "mask/{sample}_repeatmasker_{n}.fa.tbl"
  params:
    outdir = "mask"
  threads: 2
  singularity:
    "shub://tpall/repeatmasker-singularity"
  shell:
    """
    RepeatMasker -qq -pa {threads} {input.fa} -dir {params.outdir}
    if head -n 1 {output.out} | grep -q "There were no repetitive sequences detected"
      then ln -sr {input.fa} {output.masked} \
           && touch {output.tbl}
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
    masked_filt = temp("mask/{sample}_repmaskedgood_{n}.fa"),
    original_filt = temp("mask/{sample}_unmaskedgood_{n}.fa")
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/filter/masked"

# MegaBlast against reference genome to remove host sequences
rule megablast_refgenome:
    input:
      query = rules.repeatmasker_good.output.masked_filt
    output:
      out = temp("blast/{sample}_megablast_{n}.tsv")
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

# Filter megablast records for the cutoff value
rule parse_megablast:
    input:
      blast_result = rules.megablast_refgenome.output.out,
      query = rules.repeatmasker_good.output.masked_filt
    output:
      mapped = temp("blast/{sample}_refgenome_megablast_{n}_known-host.tsv"),
      unmapped = temp("blast/{sample}_refgenome_megablast_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Collect stats from preprocess outputs.
rule preprocess_stats:
  input:
    rules.preprocess.output.trimmed,
    rules.unmapped_fasta.output,
    expand("blast/{{sample}}_refgenome_megablast_{n}_unmapped.fa", n = N),
    rules.cd_hit.output.repres,
    rules.tantan.output,
    rules.tantan_good.output,
    expand(["mask/{{sample}}_repmaskedgood_{n}.fa", "mask/{{sample}}_unmaskedgood_{n}.fa"], n = N)
  output:
    "stats/{sample}_preprocess.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]

# Refgenome mapping stats.
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
