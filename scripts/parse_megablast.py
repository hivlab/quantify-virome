# This script accepts a MegaBLAST output file that were blasted against human
# genome, find out whether the best hit has a e value lower than the cutoff. If
# yes, output query information. If no, the sequence will be kept for further analysis.
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93
from Bio import SearchIO
blast_results = SearchIO.parse(snakemake.input[0], 'blast-xml')

with open(snakemake.output[0], "w") as out:
  for blast_qresult in blast_results:
    if len(blast_qresult) > 0:
        if blast_qresult[0][0].evalue < snakemake.params["e_cutoff"]:
          print(blast_qresult[0][0], file = out)
