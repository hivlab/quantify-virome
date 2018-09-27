
from Bio import SearchIO
from Bio import SeqIO
import os
import pandas as pd


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

# https://www.biostars.org/p/10162/
def get_ids(source):
    ids = map(lambda x: x.id, source)
    return set(ids)

def subset_records(source, ids):
    for record in source:
        if record.id in ids:
            yield record


def parse_blast(blast_xml, unknowns_masked_fasta, known_out_xml, unknown_out_fasta, evalue_threshold = 1e-10):
    """Parse blast output
    Write known sequences to file, return unmapped ids
    - blast_xml blast-xml file.
    - unknowns_masked_fasta repeatmasker output fasta file to be subset.
    - known_out_xml significant hits (evalue < evalue_threshold) from blast_xml in xml format.
    - unknown_out_fasta unmapped hits (evalue > evalue_threshold) from blast-xml in fasta format.
    - evalue_threshold default 1e-10
    """
    if os.stat(blast_xml).st_size == 0:
      touch(known_out_xml)
      touch(unknown_out_fasta)
      return
    blast_results = SearchIO.parse(blast_xml, 'blast-xml')
    unmapped_masked = SeqIO.index(unknowns_masked_fasta, "fasta")
    with open(known_out_xml, "w") as known_xml, open(unknown_out_fasta, "w") as unknown_fa:
        for query in blast_results:
          if query and query[0][0].evalue < evalue_threshold:
            SearchIO.write(query, known_xml, "blast-xml")
          else:
            SeqIO.write(unmapped_masked[query.id], unknown_fa, "fasta")

def subset_unmasked_csv(blast_csv, unmasked_fasta, output):
    # blast hits
    hits = pd.read_csv(blast_csv)
    # keep query id until first whitespace
    queryid_with_hits = hits["query"].str.extract("(^[^\\s]+)", expand = False)
    # convert to list
    queryid_with_hits = list(queryid_with_hits)
    # reduce to unique ids
    queryid_with_hits = set(queryid_with_hits)
    # unmasked sequences
    raw_qresults = (qresult for qresult in SeqIO.parse(unmasked_fasta, "fasta"))
    # filter unmasked using blast hits
    filtered_records = (qresult for qresult in raw_qresults if qresult.id in queryid_with_hits)
    # write to fasta
    SeqIO.write(filtered_records, output, "fasta")


######################################


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

def touch(path):
  with open(path, 'a'):
        os.utime(path, None)

# parse blast output, write known sequences to file, return unmapped ids
# blast_xml blast-xml file
# unknowns_masked_fasta repeatmasker output fasta file to be subset
# known_out_xml significant hits (evalue < evalue_threshold) from blast_xml in xml format
# unknown_out_fasta unmapped hits (evalue > evalue_threshold) from blast-xml in fasta format
# evalue_threshold default 1e-10
def parse_blast(blast_xml, unknowns_masked_fasta, known_out_xml, unknown_out_fasta, evalue_threshold = 1e-10):
  if os.stat(blast_xml).st_size == 0:
    touch(known_out_xml)
    touch(unknown_out_fasta)
    return
  blast_results = SearchIO.parse(blast_xml, 'blast-xml')
  unmapped_masked = SeqIO.index(unknowns_masked_fasta, "fasta")
  with open(known_out_xml, "w") as known_xml, open(unknown_out_fasta, "w") as unknown_fa:
      for query in blast_results:
        if query and query[0][0].evalue < evalue_threshold:
          SearchIO.write(query, known_xml, "blast-xml")
        else:
          SeqIO.write(unmapped_masked[query.id], unknown_fa, "fasta")

# subset unmasked fasta sequences using blast hits
def subset_unmasked(blast_xml, unmasked_fasta, output):
    # blast hits
    hits = SearchIO.index(blast_xml, 'blast-xml')
    # unmasked sequences
    raw_qresults = (qresult for qresult in SeqIO.parse(unmasked_fasta, "fasta"))
    # filter unmasked using blast hits
    filtered_records = (qresult for qresult in raw_qresults if qresult.id in hits)
    # write to fasta
    SeqIO.write(filtered_records, output, 'fasta')

def subset_unmasked_csv(blast_csv, unmasked_fasta, output):
    # blast hits
    hits = pd.read_csv(blast_csv)
    # keep query id until first whitespace
    queryid_with_hits = hits["query"].str.extract("(^[^\\s]+)", expand = False)
    # convert to list
    queryid_with_hits = list(queryid_with_hits)
    # reduce to unique ids
    queryid_with_hits = set(queryid_with_hits)
    # unmasked sequences
    raw_qresults = (qresult for qresult in SeqIO.parse(unmasked_fasta, "fasta"))
    # filter unmasked using blast hits
    filtered_records = (qresult for qresult in raw_qresults if qresult.id in queryid_with_hits)
    # write to fasta
    SeqIO.write(filtered_records, output, "fasta")

# subset masked fasta seqs using ref genome nonmapping fasta seq ids
def subset_masked(unmapped, masked, output):
  # unmapped
  unmapped = SeqIO.index(unmapped, "fasta")
  # masked sequences
  masked = (sequence for sequence in SeqIO.parse(masked, "fasta"))
  # filter masked using unmapped ids
  filtered_records = (sequence for sequence in masked if qresult.id in unmapped)
  # write to fasta
  SeqIO.write(filtered_records, output, 'fasta')
