
from Bio import SeqIO

def filter_N(masked, original, masked_filt, original_filt, min_length, por_n):
  """Filters out sequences with many N-s.
  Filters out sequeces with less than or equal to min_length non-N bases
  and sequences with more than por_n % N-s.
  """
  results = SeqIO.parse(masked, 'fasta')
  original = SeqIO.index(original, "fasta")
  with open(masked_filt, "w") as masked, open(original_filt, "w") as unmasked:
    for sequence in results:
      if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and
        (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
        SeqIO.write(sequence, masked, "fasta")
        SeqIO.write(original[sequence.id], unmasked, "fasta")

filter_N(snakemake.input[0], snakemake.input[1], snakemake.output[0], snakemake.output[1])
