
from Bio import SeqIO

# https://www.biostars.org/p/10162/
def get_ids(source):
    ids = map(lambda x: x.id, source)
    return set(ids)

def subset_records(source, ids):
    for record in source:
        if record.id in ids:
            yield record

# Ref genome unmapped reads
unmapped = SeqIO.parse(snakemake.input[0], "fasta")

# Get unmapped sequence ids
unmapped_ids = get_ids(unmapped)

# Subset masked sequences using unmapped_ids
masked = subset_records(SeqIO.parse(snakemake.input[1], "fasta"), unmapped_ids)

# Write unmasked subset to file
masked_count = SeqIO.write(masked, snakemake.output[0], 'fasta')
