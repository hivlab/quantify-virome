
from Bio import SeqIO
from helpers import get_ids
from helpers import subset_records

# Ref genome unmapped reads
unmapped = SeqIO.parse(snakemake.input[0], "fasta")

# Get unmapped sequence ids
unmapped_ids = get_ids(unmapped)

# Subset masked sequences using unmapped_ids
masked = subset_records(SeqIO.parse(snakemake.input[1], "fasta"), unmapped_ids)

# Write unmasked subset to file
masked_count = SeqIO.write(masked, snakemake.output[0], 'fasta')

with open(snakemake.output[0] + ".log", "w") as text_file:
  print(f"Number of sequences before filtering: {len(list(SeqIO.parse(snakemake.input[0], 'fasta')))}\nNumber of sequences after filtering: {masked_count}", file = text_file)

