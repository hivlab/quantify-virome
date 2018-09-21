
from Bio import SeqIO

def filter_N(masked, masked_filt, min_length, por_n, original = None, original_filt = None):
  """Filters out sequences with many N-s.
  Filters out sequeces with less than or equal to min_length non-N bases
  and sequences with more than por_n % N-s.
  """
  ma = SeqIO.parse(str(masked), "fasta")
  if (original is not None and original_filt is not None):
    original = SeqIO.index(str(original), "fasta")
    with open(masked_filt, "w") as maf, open(original_filt, "w") as orf:
      for record in ma:
        sequence = str(record.seq).upper()
        if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
          SeqIO.write(record, maf, "fasta")
          SeqIO.write(original[record.id], orf, "fasta")
  else:
    with open(masked_filt, "w") as maf:
      for record in ma:
        sequence = str(record.seq).upper()
        if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
          SeqIO.write(record, maf, "fasta")

def run_filter_n(input, output, params):
  # merge function arguments into dictionary
  options = dict(input)
  options.update(output)
  options.update(params)
  # run filter
  filter_N(**options)

run_filter_n(snakemake.input, snakemake.output, snakemake.params)
