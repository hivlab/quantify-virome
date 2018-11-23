
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
  urls = SAMPLES.loc[wildcards.sample, ['fq1', 'fq2']]
  return list(urls)

## Preprocessing for fastq files
# Imports local or remote fastq(.gz) files
# Downsamples runs based on user-provided fractions in samples.tsv file
# Adapter trimming
# Quality filtering
rule fastp:
  input:
    lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close = True) if config["remote"] else get_fastq
  output:
    pair1 = "munge/{sample}_pair1_trimmed.fq.gz",
    pair2 = "munge/{sample}_pair2_trimmed.fq.gz",
    html = "munge/{sample}_fastp_report.html",
    json = "munge/{sample}_fastp_report.json",
    sub1 = temp("munge/{sample}_sub1.fq.gz"),
    sub2 = temp("munge/{sample}_sub2.fq.gz")
  params:
    frac = lambda wildcards: SAMPLES.loc[wildcards.sample, ['frac']][0],
    seed = config["seed"],
    fastp = "--trim_front1 5 --trim_tail1 5 --length_required 50 --low_complexity_filter --complexity_threshold 8"
  threads: 8
  log: "logs/{sample}_fastp.log"
  conda:
      "../envs/fastp.yaml"
  shell:
    """
    if (( $(echo "{params.frac} > 0" | bc) )) && (( $(echo "{params.frac} < 1" | bc) )); then
      seqtk sample -s{params.seed} {input[0]} {params.frac} > {output.sub1}
      seqtk sample -s{params.seed} {input[1]} {params.frac} > {output.sub2}
    else
      ln -sr {input[0]} {output.sub1}
      ln -sr {input[1]} {output.sub2}
    fi
    fastp -i {output.sub1} -I {output.sub2} -o {output.pair1} -O {output.pair2} {params.fastp} -h {output.html} -j {output.json} -w {threads} > {log} 2>&1
    """

## Stitch paired reads
rule fastq_join:
  input:
    rules.fastp.output.pair1,
    rules.fastp.output.pair2
  output:
    "munge/{sample}_un1.fq.gz",
    "munge/{sample}_un2.fq.gz",
    "munge/{sample}_join.fq.gz"
  params:
    diff = config["fastq-join"]["maximum_difference"],
    overlap = config["fastq-join"]["minimum_overlap"],
    template = "munge/{sample}_%.fq.gz"
  conda:
    "../envs/fastq-join.yaml"
  log:
    "logs/{sample}_fastq_join.log"
  shell:
    """
    fastq-join -p {params.diff} -m {params.overlap} {input} -o {params.template} > {log} 2>&1
    """
