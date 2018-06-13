
## Tantan mask of low complexity DNA sequences [6]
rule tantan:
  input:
    os.path.join(config["outdir"], "{sample}/05_cdhit/merged.cdhit.fa")
  output:
    os.path.join(config["outdir"], "{sample}/06_tantan/cdhit.tantan.fa")
  params:
    "-x N"
  conda:
      "../envs/tantan.yml"
  shell:
    """
    tantan {params} {input} > {output}
    """

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
  conda:
      "../envs/biopython.yml"
  script:
      "../scripts/tantan_good.py"

## Split reads to smaller files for Repeatmasker [8]
rule split_fasta:
  input: os.path.join(config["outdir"], "{sample}/07_tantan_good/tantan.goodseq.fa")
  output:
    os.path.join(config["outdir"], dynamic("{sample}/08_split_fasta/tantan.goodseq.{n}.fa"))
  params:
    config["split_fasta"]["batch_size"],
    os.path.join(config["outdir"], "{sample}/08_split_fasta/tantan.goodseq.%i.fa")
  conda:
      "../envs/biopython.yml"
  script:
    "../scripts/split_fasta.py"

## Repeatmasker [9]
# Set RepBase library location environment variable and copy repeatmasker configuration file

shell(
"""
export REPEATMASKER_REPBASE_FILE=config["repbase_file"]
if [ ! -n "$(find $CONDA_PREFIX/share/RepeatMasker/ -maxdepth 1 -name 'RepeatMaskerConfig.pm' -print -quit)" ]
then
  cp envs/RepeatMaskerConfig.pm $CONDA_PREFIX/share/RepeatMasker/
fi
"""
)

rule repeatmasker:
  input:
    fa = os.path.join(config["outdir"], dynamic("{sample}/08_split_fasta/tantan.goodseq.{n}.fa")),
    repbase = config["repbase_file"]
  output:
    os.path.join(config["outdir"], dynamic("{sample}/09_repeatmasker/tantan.goodseq.{n}.fa.masked"))
  params:
    cluster = "-cwd -V",
    dir = os.path.join(config["outdir"], "{sample}")
  threads: 8
  shell:
    """
    RepeatMasker -qq -pa {threads} {input.fa} -dir {params.dir}/09_repeatmasker
    """

## Filter repeatmasker output [10]
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
  conda:
    "../envs/biopython.yml"
  script:
    "../scripts/repeatmasker_good.py"
