
# This script is based on https://stackoverflow.com/a/13886498/1657871

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

# Filter sequences
filtered_seq = filter_records(SeqIO.parse(snakemake.input[0], 'fasta'),
                          min_length = snakemake.params["min_length"],
                          por_n = snakemake.params["por_n"])

# Write filtered sequences to snakemake output
SeqIO.write(filtered_seq, snakemake.output[0], 'fasta')
