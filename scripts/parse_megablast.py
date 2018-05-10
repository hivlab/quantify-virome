from Bio.Blast import NCBIXML
result_handle = open(snakemake.input[0])
blast_records = NCBIXML.parse(result_handle)

with open(snakemake.output[0], "w") as out:
    for blast_record in blast_records:
       for alignment in blast_record.alignments:
             for hsp in alignment.hsps:
                if hsp.expect < snakemake.params["e_cutoff"]:
                     out.write("****Alignment****" + "\n")
                     out.write("sequence:", alignment.title + "\n")
                     out.write("length:", alignment.length + "\n")
                     out.write("e value:", hsp.expect + "\n")
                     out.write(hsp.query[0:75] + "...\n")
                     out.write(hsp.match[0:75] + "...\n")
                     out.write(hsp.sbjct[0:75] + "...\n")
