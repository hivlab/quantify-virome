
## Filter repeatmasker output
# 1) Sequences that do not have greater than 50 nt of consecutive
# sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule repeatmasker_good:
  input:
    masked = os.path.join(config["outdir"], dynamic("{sample}/09_repeatmasker/tantan.goodseq.{n}.fa.masked")),
    unmasked = os.path.join(config["outdir"], dynamic("{sample}/08_split_fasta/tantan.goodseq.{n}.fa"))
  output:
    masked = os.path.join(config["outdir"], dynamic("{sample}/10_repeatmasker_good/masked.{n}.fa")),
    unmasked = os.path.join(config["outdir"], dynamic("{sample}/10_repeatmasker_good/unmasked.{n}.fa"))
  params:
    min_length = config["repeatmasker_good"]["min_length"],
    por_n = config["repeatmasker_good"]["por_n"]
  script:
    "../scripts/repeatmasker_good.py"

## Repeatmasker [9]
# When no repeats are present, Repeatmasker will not create .masked file. We need to create one empty file 'manually'.
rule repeatmasker:
  input: os.path.join(config["outdir"], dynamic("{sample}/08_split_fasta/tantan.goodseq.{n}.fa"))
  output:
    os.path.join(config["outdir"], dynamic("{sample}/09_repeatmasker/tantan.goodseq.{n}.fa.masked"))
  params:
    cluster = "-cwd -V",
    dir = os.path.join(config["outdir"], "{sample}")
  threads:
    12
  shell:
    """
    RepeatMasker -pa {threads} {input}
    mv {params.dir}/08_split_fasta/*.fa.* {params.dir}/09_repeatmasker
    cd {params.dir}/09_repeatmasker
    if [ ! -n "$(find . -maxdepth 1 -name 'tantan.goodseq.*.fa.masked' -print -quit)" ]
    then
       touch tantan.goodseq.1.fa.masked
    fi
    """

## Split reads to smaller files for Repeatmasker [8]
rule split_fasta:
  input: os.path.join(config["outdir"], "{sample}/07_tantan_good/tantan.goodseq.fa")
  output:
    os.path.join(config["outdir"], dynamic("{sample}/08_split_fasta/tantan.goodseq.{n}.fa"))
  params:
    batch_size = 2000,
    stub = os.path.join(config["outdir"], "{sample}/08_split_fasta/tantan.goodseq.%i.fa")
  script:
    "../scripts/split_fasta.py"

## Filter tantan output [7]
# 1) Sequences that do not have greater than 50 nt of consecutive
# sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
  input:
    os.path.join(config["outdir"], "{sample}/06_tantan/cdhit.tantan.fa")
  output:
    os.path.join(config["outdir"], "{sample}/07_tantan_good/tantan.goodseq.fa")
  params:
    min_length = config["tantan_good"]["min_length"],
    por_n = config["tantan_good"]["por_n"]
  script:
      "../scripts/tantan_good.py"

## Tantan mask of low complexity DNA sequences [6]
rule tantan:
  input:
    os.path.join(config["outdir"], "{sample}/05_cdhit/merged.cdhit.fa")
  output:
    os.path.join(config["outdir"], "{sample}/06_tantan/cdhit.tantan.fa")
  params:
    "-x N"
  shell:
    """
    tantan {params} {input} > {output}
    """
