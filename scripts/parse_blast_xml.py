# Find out whether the MegaBLAST best hit has a e value lower than the cutoff. If
# yes, output query information. If no, the sequence will be kept for further analysis.
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93
from Bio import SearchIO
from Bio import SeqIO

def parse_blast(blast_results, unmapped_masked, known_out_xml, unknown_out_fasta, evalue_threshold = 1e-10):
  blast_results = SearchIO.parse(blast_results, 'blast-xml')
  unmapped_masked = SeqIO.index(unmapped_masked, "fasta")
  with open(known_out_xml, "w") as known_out, open(unknown_out_fasta, "w") as unknown_out:
    for query in blast_results:
      if query and query[0][0].evalue < evalue_threshold:
        SearchIO.write(query, known_out, "blast-xml")
      else:
        SeqIO.write(unmapped_masked[query.id], unknown_out, "fasta")

# Write known hits and unknowns to xml and fasta file, respectively
parse_blast(snakemake.input[0], snakemake.input[1], snakemake.output[0], snakemake.output[1], evalue_threshold = snakemake.params["e_cutoff"])
