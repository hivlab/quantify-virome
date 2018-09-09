
from Bio import SearchIO
from Bio.Blast import NCBIXML
from Bio import SeqIO
from helpers import subset_unmasked

subset_unmasked(blast_xml = snakemake.input[0],
           unmasked_fasta = snakemake.input[2],
                   output = snakemake.output[0])

subset_unmasked(blast_xml = snakemake.input[1],
           unmasked_fasta = snakemake.input[2],
                   output = snakemake.output[1])
