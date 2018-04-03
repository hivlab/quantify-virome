# Function to filter sequences

# min_length =  Minimum sequence length after Ns
# por_n = Maximum percent of masked bases (Ns)

def filter_records(source, min_length, por_n):
    for seq_record in source:
        sequence = str(seq_record.seq).upper()
        if ((float(len(sequence)) - float(sequence.count("N"))) >= min_length and
      (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
             yield seq_record
