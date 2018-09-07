from Bio.Blast.Applications import NcbiblastxCommandline

cline = NcbiblastxCommandline(query = snakemake.input[1],
                              db = snakemake.input[0],
                              show_gis = True,
                              num_threads = snakemake.threads,
                              evalue = snakemake.params["evalue"],
                              db_soft_mask = snakemake.params["db_soft_mask"],
                              outfmt = 5,
                              out = snakemake.output)
stdout, stderr = cline()
