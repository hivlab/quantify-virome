
from Bio import SeqIO
from helpers import batch_iterator

record_iter = SeqIO.parse(snakemake.input[0], "fasta")
for i, batch in enumerate(batch_iterator(record_iter, snakemake.params["batch_size"])):
    filename = snakemake.params["stub"] % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
