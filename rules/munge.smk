
def get_fastq(wildcards):
 fastq_files = samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
 paths = []
 for x in fastq_files:
   path = os.path.join(config["datadir"], x)
   paths.append(path)
 return paths

## Run cd-hit to find and munge duplicate reads [5]
rule cd_hit:
  input:
    os.path.join(config["outdir"], "{sample}/04_fasta/stitched.merged.fasta")
  output:
    clusters = os.path.join(config["outdir"], "{sample}/05_cdhit/merged.cdhit.fa"),
    report = os.path.join(config["outdir"], "{sample}/05_cdhit/stitched.merged.cdhit.report")
  params:
    "-c 0.984 -G 0 -n 8 -d 0 -aS 0.984 -g 1 -r 1 -M 0"
  threads:
    config["cd-hit"]["threads"]
  conda:
    "../envs/cd-hit.yml"
  shell:
    """
    cd-hit-est -i {input} -o {output.clusters} {params} -T {threads} {params} > {output.report}
    """

## Convert fastq to fasta format [4]
rule fastq2fasta:
  input: os.path.join(config["outdir"], "{sample}/03_merged/stitched.merged.fq.gz")
  output: os.path.join(config["outdir"], "{sample}/04_fasta/stitched.merged.fasta")
  shell:
    """
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
    """

## Merge stitched reads
rule merge_reads:
  input:
    lambda wildcards: expand(os.path.join(config["outdir"], "{sample}/02_stitched/{pair}.fq.gz"), sample = wildcards.sample, pair = ["join", "un1", "un2"])
  output:
    os.path.join(config["outdir"], "{sample}/03_merged/stitched.merged.fq.gz")
  shell:
    """
    cat {input[0]} {input[1]} {input[2]} > {output[0]}
    """

## Stitch paired reads [2]
rule fastq_join:
  input:
    lambda wildcards: expand(os.path.join(config["outdir"], "{sample}/01_fastp/pair{n}.truncated.gz"), sample = wildcards.sample, n = [1, 2])
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

## All-in-one preprocessing for FastQ files [1,3]
# Adapter trimming is enabled by default
# Quality filtering is enabled by default
# Replaces AdapteRemoval, prinseq and fastqc
rule fastp:
    input:
      get_fastq
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
      fastp -i {input[0]} -I {input[1]} -o {output.pair1} -O {output.pair2} {params.options} -h {params.html} -j {params.json} -w {threads}
      """
