# Find out whether the MegaBLAST best hit has a e value lower than the cutoff. If
# yes, output query information. If no, the sequence will be kept for further analysis.
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93
from Bio import SearchIO
from Bio import SeqIO

def parse_blast(blast_xml, unknowns_masked_fasta, known_out_xml, known_out_fasta, unknown_out_fasta, evalue_threshold = 1e-10, unknowns_fasta = None):
  blast_results = SearchIO.parse(blast_xml, 'blast-xml')
  unmapped_masked = SeqIO.index(unknowns_masked_fasta, "fasta")
  if unknowns_fasta is not None:
    unmapped = SeqIO.index(unknowns_fasta, "fasta")
  with open(known_out_xml, "w") as known_xml, open(known_out_fasta, "w") as known_fa, open(unknown_out_fasta, "w") as unknown_fa:
    for query in blast_results:
      if query and query[0][0].evalue < evalue_threshold:
        SearchIO.write(query, known_xml, "blast-xml")
        if unknowns_fasta is not None:
          SeqIO.write(unmapped[query.id], known_fa, "fasta")
      else:
        SeqIO.write(unmapped_masked[query.id], unknown_fa, "fasta")

# Write known hits and unknowns to xml and fasta file, respectively
parse_blast(blast_xml = snakemake.input[0],
            unknowns_masked_fasta = snakemake.input[1],
            known_out_xml = snakemake.output[0],
          known_out_fasta = snakemake.output[1],
        unknown_out_fasta = snakemake.output[2],
        evalue_threshold = snakemake.params["e_cutoff"])
