
from Bio import SeqIO
import pandas as pd

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

subset_unmasked_csv(blast_csv = snakemake.input[0],
           unmasked_fasta = snakemake.input[1],
                   output = snakemake.output[0])
