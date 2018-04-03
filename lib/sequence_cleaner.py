
# This script is based on https://stackoverflow.com/a/13886498/1657871

from Bio import SeqIO
from lib import filter_records

# Filter sequences
input_seq = SeqIO.parse(snakemake.input[0], 'fasta')
filtered_seq = filter_records(input_seq,
                          min_length = snakemake.params["min_length"],
                          por_n = snakemake.params["por_n"])

# Write filtered sequences to snakemake output
SeqIO.write(filtered_seq, snakemake.output[0], 'fasta')
