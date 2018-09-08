
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
# blast_xml blast-xml file
# unknowns_masked_fasta repeatmasker output fasta file to be subset
# known_out_xml significant hits (evalue < evalue_threshold) from blast_xml in xml format
# unknown_out_fasta unmapped hits (evalue > evalue_threshold) from blast-xml in fasta format
# evalue_threshold default 1e-10
# unknowns_fasta unmasked fasta file to be subset, default None
# known_out_fasta significant hits (evalue < evalue_threshold) from blast_xml in fasta format, default None
def parse_blast(blast_xml, unknowns_masked_fasta, known_out_xml, unknown_out_fasta, evalue_threshold = 1e-10, unknowns_fasta = None, known_out_fasta = None):
  blast_results = SearchIO.parse(blast_xml, 'blast-xml')
  unmapped_masked = SeqIO.index(unknowns_masked_fasta, "fasta")
  if unknowns_fasta is not None:
    if known_out_fasta is None: raise Exception("Please provide known_out_fasta file path.")
    unmapped = SeqIO.index(unknowns_fasta, "fasta")
    with open(known_out_xml, "w") as known_xml, open(known_out_fasta, "w") as known_fa, open(unknown_out_fasta, "w") as unknown_fa:
      for query in blast_results:
        if query and query[0][0].evalue < evalue_threshold:
          SearchIO.write(query, known_xml, "blast-xml")
          SeqIO.write(unmapped[query.id], known_fa, "fasta")
        else:
          SeqIO.write(unmapped_masked[query.id], unknown_fa, "fasta")
  else:
    with open(known_out_xml, "w") as known_xml, open(unknown_out_fasta, "w") as unknown_fa:
      for query in blast_results:
        if query and query[0][0].evalue < evalue_threshold:
          SearchIO.write(query, known_xml, "blast-xml")
        else:
          SeqIO.write(unmapped_masked[query.id], unknown_fa, "fasta")
