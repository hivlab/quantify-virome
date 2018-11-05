
from Bio import SeqIO

def filter_N(masked, masked_filt, min_length, por_n, fa = None, fa_filt = None):
  """Filters out sequences with many N-s.
  Filters out sequeces with less than or equal to min_length non-N bases
  and sequences with more than por_n % N-s.
  """
  ma = SeqIO.parse(str(masked), "fasta")
  if (fa is not None and fa_filt is not None):
    fa = SeqIO.index(str(fa), "fasta")
    with open(masked_filt, "w") as maf, open(fa_filt, "w") as orf:
      for record in ma:
        sequence = str(record.seq).upper()
        if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
          SeqIO.write(record, maf, "fasta")
          SeqIO.write(fa[record.id], orf, "fasta")
  else:
    with open(masked_filt, "w") as maf:
      for record in ma:
        sequence = str(record.seq).upper()
        if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
          SeqIO.write(record, maf, "fasta")

def run_filter_n(input, output, params):
  # Merge function arguments into dictionary
  options = dict(input)
  options.update(output)
  options.update(params)

  # Keep only filter_N arguments
  args = {'masked', 'masked_filt', 'min_length', 'por_n', 'fa', 'fa_filt'}
  options = {k: options[k] for k in options.keys() & args}

  # Run filter
  filter_N(**options)
