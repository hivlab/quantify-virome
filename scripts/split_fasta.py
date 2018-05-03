
from Bio import SeqIO
from helpers import batch_iterator

record_iter = SeqIO.parse(snakemake.input[0], "fasta")
for i, batch in enumerate(batch_iterator(record_iter, snakemake.params[0])):
    filename = snakemake.params[1] % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
