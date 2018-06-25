# Find out whether the MegaBLAST best hit has a e value lower than the cutoff. If
# yes, output query information. If no, the sequence will be kept for further analysis.
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93
from Bio import SearchIO
from Bio import SeqIO
from helpers import parse_blast
from helpers import subset_records

# Write known hits to file and return unmapped ids
keep = parse_blast(blastxml = snakemake.input["blastxml"], evalue_threshold = snakemake.params["e_cutoff"], outfile = snakemake.output["known"])
keep_lst = list(keep)

# Subset masked sequences using unmapped_ids
unmapped_masked = SeqIO.parse(snakemake.input["query"], "fasta")
unmapped_keep = subset_records(unmapped_masked, keep_lst)

# Write unmasked subset to file
masked_count = SeqIO.write(unmapped_keep, snakemake.output["unmapped"], "fasta")
