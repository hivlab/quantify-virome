
def filter_records(source, min_length, por_n):
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

# http://biopython.org/wiki/Split_large_file

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

# parse blast output, write known sequences to file, return unmapped ids
def parse_blast(blast_results, outfile, evalue_threshold = 1e-10):
  with open(outfile, "w") as known:
    for query in blast_results:
      if query and query[0][0].evalue < evalue_threshold):
        known.write("%s\n\n" % query[0][0])
      else:
        yield query.id
