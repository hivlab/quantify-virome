# Find out whether the MegaBLAST best hit has a e value lower than the cutoff. If
# yes, output query information. If no, the sequence will be kept for further analysis.
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93
from helpers import parse_blast

# Write known hits and unknowns to xml and fasta file, respectively
parse_blast(blast_xml = snakemake.input[0], unknowns_masked_fasta = snakemake.input[1], known_out_xml = snakemake.output[0], unknown_out_fasta = snakemake.output[1], evalue_threshold = snakemake.params["e_cutoff"])
