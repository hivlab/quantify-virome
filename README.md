
[![Travis-CI Build Status](https://travis-ci.org/<USERNAME>/<REPO>.svg?branch=master)](https://travis-ci.org/<USERNAME>/<REPO>)

The goal of this repo is to reproducibly recreate VirusSeeker Virome workflow.

# VirusSeeker Virome workflow

Following workflow description is a copy-paste from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5326578/

## Sequence preprocessing
The preprocessing of sequence files consists of the following steps: 

1. **Trim adapter and/or primer sequences.**

2. **Join read1 and 2 (of paired end reads) together to form a longer read** if they overlap by defined criteria using fastq-join in the ea-utils package (https://github.com/ExpressionAnalysis/ea-utils) [55].

3. **Quality filtering of reads** (trim low quality nucleotides, poly A sequences, remove reads with low average quality score) to obtain good quality sequences using PRINSEQ [56].

4. **Remove redundant sequences.** Identical or nearly-identical sequences are frequently present in NGS data, either due to the sheer depth of NGS or because many of the pre-sequencing sample preparation methods involve PCR amplification. To reduce the computing cost of downstream analysis and to reduce amplification artifacts, CD-HIT [57] is used to cluster similar sequences. **The default parameters in VS-Virome are set to cluster sequences that share ≥ 95% identity over 95% of the sequence length.** The longest sequence from each cluster is retained as the representative sequence and used for downstream analysis. These are the “unique sequences”. 

5. **Mask repetitive sequences and sequence quality control.** Many eukaryotic genomes contain stretches of highly repetitive DNA sequences which cause problems in BLAST-based similarity searches and result in high rates of false-positive alignments. Tantan [58] and RepeatMasker (http://www.repeatmasker.org) [59] are used to mask interspersed repeats and low complexity DNA sequences. A sequence fails the quality control criteria if it does not contain a stretch of at least 50 consecutive non-“N” nucleotides (i.e., “Filtered sequence”) or if greater than 40% of the total length of the sequence is masked (i.e., “low complexity sequence”). These sequences are removed from further analysis. Remaining sequences are “good sequences”.

6. **Remove host sequences by aligning sequences to reference genome using BWA-MEM [60] and MegaBLAST.** Any sequence mapped to “Host” genomic sequence is removed from further analysis. After sequence preprocessing we obtain high quality, unique, non-host sequences.


## Workflow graph

![Virome workflow](graph/dag.svg)

## Installation


## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```
