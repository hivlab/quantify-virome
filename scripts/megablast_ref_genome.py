from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(query = snakemake.input[1],
                                db = snakemake.input[0],
                                num_threads = snakemake.threads,
                                perc_identity = snakemake.params["perc_ident"],
                                evalue = snakemake.params["evalue"],
                                word_size = snakemake.params["word_size"],
                                num_descriptions = snakemake.params["num_desc"],
                                num_alignments = snakemake.params["num_align"],
                                outfmt = 5,
                                out = snakemake.output)
stdout, stderr = blastn_cline()
