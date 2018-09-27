
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline

def run_blast(input, output, params):
  # merge function arguments into dictionary
  options = dict(input)
  options.update(output)
  options.update(params)
  # compose blastn or blastx command
  if "task" in options:
    cline = NcbiblastnCommandline(**options)
  else:
    cline = NcbiblastxCommandline(**options)
  # run blast
  stdout, stderr = cline()

run_blast(snakemake.input, snakemake.output, snakemake.params)
