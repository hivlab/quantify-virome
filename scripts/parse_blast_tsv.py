# Find out whether the MegaBLAST best hit has a e value lower than the cutoff. If
# yes, output query information. If no, the sequence will be kept for further analysis.
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93
from common.helpers import parse_blast_fmt6

def run_parse_blast(input, output, params):
  # merge function arguments into dictionary
  options = dict(input)
  options.update(output)
  options.update(params)
  # unwrap arguments and run function
  parse_blast_fmt6(**options)

run_parse_blast(snakemake.input, snakemake.output, snakemake.params)
