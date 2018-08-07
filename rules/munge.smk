
def get_fastq(wildcards, path = '.', read_pair='fq1'):
 fq = samples.loc[wildcards.sample, [read_pair]].dropna()[0]
 return os.path.join(path, fq)

## All-in-one preprocessing for FastQ files [1,3]
# Adapter trimming is enabled by default
# Quality filtering is enabled by default
# Replaces AdapteRemoval, prinseq and fastqc
rule fastp:
    input:
      fq1 = lambda wildcards: get_fastq(wildcards, config["datadir"], 'fq1'),
      fq2 = lambda wildcards: get_fastq(wildcards, config["datadir"], 'fq2')
    output:
      pair1 = os.path.join(config["outdir"], "{sample}/01_fastp/pair1.truncated.gz"),
      pair2 = os.path.join(config["outdir"], "{sample}/01_fastp/pair2.truncated.gz")
    params:
      options = "-f 5 -t 5 -l 50 -y -Y 8",
      html = os.path.join(config["outdir"], "{sample}/01_fastp/report.html"),
      json = os.path.join(config["outdir"], "{sample}/01_fastp/report.json")
    threads:
      config["fastp"]["threads"]
    conda:
      "../envs/fastp.yml"
    shell:
      """
      fastp -i {input.fq1} -I {input.fq2} -o {output.pair1} -O {output.pair2} {params.options} -h {params.html} -j {params.json} -w {threads}
      """

## Stitch paired reads [2]
rule fastq_join:
  input: rules.fastp.output
  output:
    os.path.join(config["outdir"], "{sample}/02_stitched/join.fq.gz"),
    os.path.join(config["outdir"], "{sample}/02_stitched/un1.fq.gz"),
    os.path.join(config["outdir"], "{sample}/02_stitched/un2.fq.gz")
  params:
    config["fastq-join"]["maximum_difference"],
    config["fastq-join"]["minimum_overlap"],
    template = os.path.join(config["outdir"], "{sample}/02_stitched/%.fq.gz")
  conda:
    "../envs/fastq-join.yml"
  shell:
    """
    fastq-join \
    -p {params[0]} \
    -m {params[1]} \
    {input[0]} \
    {input[1]} \
    -o {params[2]}
    """

## Merge stitched reads [3]
rule merge_reads:
  input: rules.fastq_join.output
  output:
    os.path.join(config["outdir"], "{sample}/03_merged/stitched.merged.fq.gz")
  shell:
    """
    cat {input[0]} {input[1]} {input[2]} > {output[0]}
    """

## Convert fastq to fasta format [4]
rule fastq2fasta:
  input: rules.merge_reads.output
  output: os.path.join(config["outdir"], "{sample}/04_fasta/stitched.merged.fasta")
  shell:
    """
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
    """

## Run cd-hit to find and munge duplicate reads [5]
rule cd_hit:
  input: rules.fastq2fasta.output
  output:
    clusters = os.path.join(config["outdir"], "{sample}/05_cdhit/merged.cdhit.fa"),
    report = os.path.join(config["outdir"], "{sample}/05_cdhit/stitched.merged.cdhit.report")
  params:
    "-c 0.984 -G 0 -n 8 -d 0 -aS 0.984 -g 1 -r 1 -M 0 -T 0"
  conda:
    "../envs/cd-hit.yml"
  shell:
    """
    cd-hit-est -i {input} -o {output.clusters} {params} > {output.report}
    """
