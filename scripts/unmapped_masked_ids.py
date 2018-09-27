
from Bio import SeqIO
from common.helpers import get_ids
from common.helpers import subset_records

# Ref genome unmapped reads
unmapped = SeqIO.parse(snakemake.input[0], "fasta")

# Get unmapped sequence ids
unmapped_ids = get_ids(unmapped)

# Subset masked sequences using unmapped_ids
masked = subset_records(SeqIO.parse(snakemake.input[1], "fasta"), unmapped_ids)

# Write unmasked subset to file
masked_count = SeqIO.write(masked, snakemake.output[0], 'fasta')
