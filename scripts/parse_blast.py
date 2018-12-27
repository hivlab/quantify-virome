
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import pandas as pd
from Bio import SeqIO
from pandas.io.common import EmptyDataError
from common.helpers import touch
# Find out whether the BLAST best hit has a e value lower than the cutoff. If
# yes, output query information. If no, the sequence will be kept for further analysis.
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93
# Function expexts BLAST tabular format (outfmt 6)
def read_data(file):
    try:
        df = pd.read_table(file)
    except EmptyDataError:
        df = pd.DataFrame()
    return df

def parse_blast_result(blast_result, query, e_cutoff, outfmt, mapped, unmapped):
  # Import blast output table
  tab = read_data(blast_result)
  if len(tab.index) == 0:
    known_ids = set()
    touch(mapped)
  else:
    # Import column names, replace std when present, remove double quote
    std = "'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    # Munge outfmt string
    outfmt = outfmt.replace("std", std)
    outfmt = outfmt.replace("'", "")
    # Get colnames from outfmt string
    colnames = [colname for colname in outfmt.split() if "6" not in colname]
    # Assign column names
    tab.columns = colnames
    # Filter results
    known = tab[(tab.evalue <= e_cutoff)]
    # Write seqs below threshold to file
    known.to_csv(mapped, sep = '\t', encoding = "utf-8", index = False)
    known_ids = set(known.qseqid)
  # Subset blast input
  with open(unmapped, "w") as out:
    for record in SeqIO.parse(str(query), "fasta"):
        if record.id not in known_ids:
            SeqIO.write(record, out, "fasta")

def run_parse_blast(input, output, params):
  # Merge function arguments into dictionary
  options = dict(input)
  options.update(output)
  options.update(params)
  # Unwrap arguments and run function
  parse_blast_result(**options)

if __name__ == '__main__':
    run_parse_blast(snakemake.input, snakemake.output, snakemake.params)
