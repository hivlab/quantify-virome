
from Bio import SearchIO
from Bio import SeqIO
import os
import pandas as pd
from pandas.io.common import EmptyDataError

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

# https://stackoverflow.com/questions/1158076/implement-touch-using-python#1160227
def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

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

def read_data(file):
    try:
        df = pd.read_table(file)
    except EmptyDataError:
        df = pd.DataFrame()
    return df 

def parse_blast_fmt6(blast_result, query, e_cutoff, outfmt, known_host, unmapped):
  # import blast output table
  tab = read_data(blast_result)
  if len(tab.index) == 0:
    known_ids = set()
    touch(known_host)
  else:
    # import column names
    colnames = list(filter(lambda x: '6' not in x, outfmt.split()))
    # assign column names
    tab.columns = colnames
    # filter results
    known = tab[(tab.evalue <= e_cutoff)]
    # write seqs below threshold to file
    known.to_csv(known_host, sep = '\t', encoding = 'utf-8', index = False)
    known_ids = set(known.qseqid)
  # subset blast input
  with open(unmapped, "w") as out:
    for record in SeqIO.parse(str(query), "fasta"):
        if record.id not in known_ids:
            SeqIO.write(record, out, "fasta")

