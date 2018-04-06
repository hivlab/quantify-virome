
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
    mv *.fa.* ../repeatmasker
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
