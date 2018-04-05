
from Bio import SeqIO
from lib.helpers import filter_records, get_ids, subset_records

masked_in = "output/repeatmasker/I1164_12629_Harvard_SIV_196_06_2_24_12_mini.tantan.goodseq.3.fa.masked"
tantan_in = "output/split_fasta/I1164_12629_Harvard_SIV_196_06_2_24_12_mini.tantan.goodseq.3.fa"

masked_filt_out = "output/repeatmasker_good/I1164_12629_Harvard_SIV_196_06_2_24_12_mini.repeatmasker.goodseq.masked.3.fa"
tantan_filt_out = "output/repeatmasker_good/I1164_12629_Harvard_SIV_196_06_2_24_12_mini.repeatmasker.goodseq.unmasked.3.fa"

# Filter masked sequences
masked = SeqIO.parse(masked_in, "fasta")
masked_filt = list(filter_records(masked, min_length = 50, por_n = 40))

# Write filtered sequences to file
masked_filt_count = SeqIO.write(masked_filt, masked_filt_out, 'fasta')

# Get filtered sequence ids
masked_filt_ids = get_ids(masked_filt)

# Subset unmasked sequences using ids
tantan_filt = subset_records(SeqIO.parse(tantan_in, "fasta"), masked_filt_ids)

# Write unmasked subset to file
tantan_filt_count = SeqIO.write(tantan_filt, tantan_filt_out, 'fasta')

with open(masked_filt_out + ".log", "w") as text_file:
    print(f"Number of sequences before filtering: {len(list(SeqIO.parse(masked_in, 'fasta')))}\nNumber of sequences after filtering: {masked_filt_count}", file = text_file)

