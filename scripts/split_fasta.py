
from re import search
from math import ceil
from Bio import SeqIO
from common.helpers import batch_iterator

# Number of sequences in fasta file
fasta_file = open(snakemake.input[0])
seqs = sum([bool(search(r"^>", line)) for line in fasta_file])

# Calculate batch size given number of files
batch_size = ceil(seqs / snakemake.params[0])

# Split sequences into chunks based on batch size and write into files
record_iter = SeqIO.parse(snakemake.input[0], "fasta")
for n, batch in enumerate(batch_iterator(record_iter, batch_size), start = 1):
    if n == snakemake.params[1]:
        with open(snakemake.output[0], "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
            break
