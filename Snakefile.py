
## Use os.path to update file paths from config file
import os.path

## Load configuration file with sample and path info
configfile: "config.yaml"

## Target rule
rule all:
    input:
      expand(os.path.join(config["outdir"], "repeatmasker_good/{sample}.repeatmasker.goodseq.{mask}.{n}.fa"), sample = config["samples"], n = [1,2,3], mask = ["masked", "unmasked"])

## Filter repeatmasker output
# 1) Sequences that do not have greater than 50 nt of consecutive
# sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule repeatmasker_good:
  input:
    masked = os.path.join(config["outdir"], dynamic("repeatmasker/{sample}.tantan.goodseq.{n}.fa.masked")),
    unmasked = os.path.join(config["outdir"], dynamic("split_fasta/{sample}.tantan.goodseq.{n}.fa"))
  output:
    masked = os.path.join(config["outdir"], dynamic("repeatmasker_good/{sample}.repeatmasker.goodseq.masked.{n}.fa")),
    unmasked = os.path.join(config["outdir"], dynamic("repeatmasker_good/{sample}.repeatmasker.goodseq.unmasked.{n}.fa"))
  params:
    min_length = 50,
    por_n = 40
  script:
    "scripts/repeatmasker_good.py"

## Repeatmasker [9]
rule repeatmasker:
  input: os.path.join(config["outdir"], dynamic("split_fasta/{sample}.tantan.goodseq.{n}.fa"))
  output:
    os.path.join(config["outdir"], dynamic("repeatmasker/{sample}.tantan.goodseq.{n}.fa.masked"))
  params:
    cluster = "-cwd -V",
    dir = os.path.join(config["outdir"], "split_fasta")
  threads:
    12
  shell:
    """
    RepeatMasker -pa {threads} {input}
    cd {params.dir}
    mv *.fa.* ../repeatmasker/
    """

## Split reads to smaller files for Repeatmasker [8]
rule split_fasta:
  input: os.path.join(config["outdir"], "tantan_good/{sample}.tantan.goodseq.fa")
  output:
    os.path.join(config["outdir"], dynamic("split_fasta/{sample}.tantan.goodseq.{n}.fa"))
  params:
    batch_size = 2000,
    stub = os.path.join(config["outdir"], "split_fasta/{sample}.tantan.goodseq.%i.fa")
  script:
    "scripts/split_fasta.py"

## Filter tantan output [7]
# 1) Sequences that do not have greater than 50 nt of consecutive
# sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
  input:
    os.path.join(config["outdir"], "tantan/{sample}.stitched.merged.cdhit.tantan.fa")
  output:
    os.path.join(config["outdir"], "tantan_good/{sample}.tantan.goodseq.fa")
  params:
    min_length = 50,
    por_n = 40
  script:
      "scripts/tantan_good.py"

## Tantan mask of low complexity DNA sequences [6]
rule tantan:
  input:
    os.path.join(config["outdir"], "cdhit/{sample}.stitched.merged.cdhit.fa")
  output:
    os.path.join(config["outdir"], "tantan/{sample}.stitched.merged.cdhit.tantan.fa")
  params:
    "-x N"
  shell:
    """
    tantan {params} {input} > {output}
    """

## Run cd-hit to find and munge duplicate reads [5]
rule cd_hit:
  input:
    os.path.join(config["outdir"], "fasta/{sample}.stitched.merged.fasta")
  output:
    clusters = os.path.join(config["outdir"], "cdhit/{sample}.stitched.merged.cdhit.fa"),
    report = os.path.join(config["outdir"], "cdhit/{sample}.stitched.merged.cdhit.report")
  resources:
    mem = 20480
  params:
    "-c 0.984 -G 0 -n 8 -d 0 -aS 0.984 -g 1 -r 1"
  threads:
    40
  shell:
    """
    cd-hit-est -i {input} -o {output.clusters} {params} -T {threads} -M {resources.mem} {params} > {output.report}
    """

## Convert fastq to fasta format [4]
rule fastq2fasta:
  input: os.path.join(config["outdir"], "merged/{sample}.stitched.merged.fq.gz")
  output: os.path.join(config["outdir"], "fasta/{sample}.stitched.merged.fasta")
  shell:
    """
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
    """


## Merge stitched reads
rule merge_reads:
  input:
    join = os.path.join(config["outdir"], "stitched/{sample}.join.fq.gz"),
    un1 = os.path.join(config["outdir"], "stitched/{sample}.un1.fq.gz"),
    un2 = os.path.join(config["outdir"], "stitched/{sample}.un2.fq.gz")
  output:
    merged = os.path.join(config["outdir"], "merged/{sample}.stitched.merged.fq.gz")
  shell:
    """
    cat {input.join} {input.un1} {input.un2} > {output.merged}
    """

## Stitch paired reads [2]
rule fastq_join:
  input:
    pair1 = os.path.join(config["outdir"], "fastp/{sample}.pair1.truncated.gz"),
    pair2 = os.path.join(config["outdir"], "fastp/{sample}.pair2.truncated.gz")
  output:
    os.path.join(config["outdir"], "stitched/{sample}.join.fq.gz"),
    os.path.join(config["outdir"], "stitched/{sample}.un1.fq.gz"),
    os.path.join(config["outdir"], "stitched/{sample}.un2.fq.gz")
  params:
    maximum_difference = 5,
    minimum_overlap = 10,
    template = os.path.join(config["outdir"], "stitched/{sample}.%.fq.gz")
  shell:
    """
    fastq-join \
    -p {params.maximum_difference} \
    -m {params.minimum_overlap} \
    {input.pair1} \
    {input.pair2} \
    -o {params.template}
    """

## All-in-one preprocessing for FastQ files [1,3]
# Adapter trimming is enabled by default
# Quality filtering is enabled by default
# Replaces AdapteRemoval, prinseq and fastqc
rule fastp:
    input:
        R1 = os.path.join(config["datadir"], "{sample}_SE1.fastq.gz"),
        R2 = os.path.join(config["datadir"], "{sample}_SE2.fastq.gz")
    output:
        pair1 = os.path.join(config["outdir"], "fastp/{sample}.pair1.truncated.gz"),
        pair2 = os.path.join(config["outdir"], "fastp/{sample}.pair2.truncated.gz")
    params:
        options = "-f 5 -t 5 -l 50 -y -Y 8",
        html = os.path.join(config["outdir"], "fastp/{sample}.report.html"),
        json = os.path.join(config["outdir"], "fastp/{sample}.report.json")
    threads: 16
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.pair1} -O {output.pair2} {params.options} -h {params.html} -j {params.json} -w {threads}
        """
