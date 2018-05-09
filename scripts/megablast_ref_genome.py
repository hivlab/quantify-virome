from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(query = input.query,
                                db = input.db,
                                num_threads = threads,
                                perc_identity = params.perc_ident,
                                evalue = params.evalue,
                                word_size = params.word_size,
                                num_descriptions = params.num_desc,
                                num_alignments = params.num_align,
                                outfmt = 5,
                                out = output)
stdout, stderr = blastn_cline()
