
from Bio import SeqIO

def filter_records(source, min_length, por_n):
  """Function to filter sequences

  min_length =  Minimum sequence length after Ns
  por_n = Maximum percent of masked bases (Ns)
  """
  for seq_record in source:
    sequence = str(seq_record.seq).upper()
    if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and
      (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
        yield seq_record

# https://www.biostars.org/p/10162/
def get_ids(source):
  ids = map(lambda x: x.id, source)
  return set(ids)

def subset_records(source, ids):
  for record in source:
    if record.id in ids:
      yield record

# Filter masked sequences
masked_filt = list(filter_records(SeqIO.parse(snakemake.input[0], "fasta"), min_length = snakemake.params["min_length"], por_n = snakemake.params["por_n"]))

# Write filtered sequences to file
masked_filt_count = SeqIO.write(masked_filt, snakemake.output[0], 'fasta')

# Get filtered sequence ids
masked_filt_ids = get_ids(masked_filt)

# Subset unmasked sequences using ids
tantan_filt = subset_records(SeqIO.parse(snakemake.input[1], "fasta"), masked_filt_ids)

# Write unmasked subset to file
tantan_filt_count = SeqIO.write(tantan_filt, snakemake.output[1], 'fasta')

with open(snakemake.output[0] + ".log", "w") as text_file:
  print(f"Number of sequences before filtering: {len(list(SeqIO.parse(snakemake.input[0], 'fasta')))}\nNumber of sequences after filtering: {masked_filt_count}", file = text_file)

