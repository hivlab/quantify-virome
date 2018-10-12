
## Tantan mask of low complexity DNA sequences [6]
rule tantan:
  input: rules.cd_hit.output.clusters
  output:
    "mask/{sample}_tantan.fasta"
  conda:
      "../envs/tantan.yml"
  shell:
    """
    tantan -x N {input} > {output}
    """

## Filter tantan output [7]
# 1) Sequences that do not have greater than 50 nt of consecutive
# sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
  input:
    masked = rules.tantan.output
  output:
    masked_filt = "mask/{sample}_tantangood.fasta"
  params:
    min_length = 50,
    por_n = 40
  conda:
      "../envs/biopython.yml"
  script:
      "../scripts/filter_masked.py"

## Split reads to smaller chunks for Repeatmasker [8]
rule split_fasta:
  input:
    rules.tantan_good.output
  output:
    "mask/{sample}_repeatmasker_{n}.fa"
  params:
    n_files,
    lambda wildcards: wildcards.n
  conda:
    "../envs/biopython.yml"
  script:
    "../scripts/split_fasta.py"

## Repeatmasker [9]
# Set RepBase library location environment variable and copy repeatmasker configuration file

shell(
"""
if [ ! -n "$(find $CONDA_PREFIX/share/RepeatMasker/ -maxdepth 1 -name 'RepeatMaskerConfig.pm' -print -quit)" ]
then
  cp envs/RepeatMaskerConfig.pm $CONDA_PREFIX/share/RepeatMasker/
fi
"""
)

# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected copy all input files to output
rule repeatmasker:
  input:
    fa = rules.split_fasta.output,
    repbase = config["repbase_file"]
  output:
    masked = "mask/{sample}_repeatmasker_{n}.fa.masked",
    out = "mask/{sample}_repeatmasker_{n}.fa.out",
    tbl = "mask/{sample}_repeatmasker_{n}.fa.tbl"
  threads: 8
  shell:
    """
    export REPEATMASKER_REPBASE_FILE={input.repbase}
    RepeatMasker -qq -pa {threads} {input.fa} -dir mask
    if head -n 1 {output.out} | grep -q "There were no repetitive sequences detected"
      then cp {input.fa} {output.masked}
    fi
    """

## Filter repeatmasker output [10]
# 1) Sequences that do not have greater than 50 nt of consecutive
# sequence without N
# 2) Sequences with >= 40% of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
  input:
    masked = rules.repeatmasker.output.masked,
    original = rules.split_fasta.output
  output:
    masked_filt = "mask/{sample}_repmaskedgood_{n}.fa",
    original_filt = "mask/{sample}_unmaskedgood_{n}.fa"
  params:
    min_length = 50,
    por_n = 40
  conda:
    "../envs/biopython.yml"
  script:
    "../scripts/filter_masked.py"
