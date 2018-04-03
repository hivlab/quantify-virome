
from Bio import SeqIO
from lib import sequence_cleaner
import os

masked = SeqIO.parse("output/repeatmasker/I1164_12629_Harvard_SIV_196_06_2_24_12_mini.tantan.goodseq.1.fa.masked", "fasta")
tantan_goodseq = SeqIO.parse("output/split_fasta/I1164_12629_Harvard_SIV_196_06_2_24_12_mini.tantan.goodseq.1.fa", "fasta")
for record in seq:
    print(record.id)



