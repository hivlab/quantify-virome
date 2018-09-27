
from common.helpers import subset_unmasked_csv

subset_unmasked_csv(blast_csv = snakemake.input[0],
           unmasked_fasta = snakemake.input[1],
                   output = snakemake.output[0])

