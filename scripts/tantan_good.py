
# This script is based on https://stackoverflow.com/a/13886498/1657871

from Bio import SeqIO
from helpers import filter_records

# Filter sequences
filtered_seq = filter_records(SeqIO.parse(snakemake.input[0], 'fasta'),
                          min_length = snakemake.params["min_length"],
                          por_n = snakemake.params["por_n"])

# Write filtered sequences to snakemake output
SeqIO.write(filtered_seq, snakemake.output[0], 'fasta')
