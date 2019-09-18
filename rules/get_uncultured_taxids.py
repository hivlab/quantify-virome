from Bio import Entrez

Entrez.email = {snakemake.params["email"]}
resp = Entrez.esearch(db = "taxonomy", term = '"environmental samples"[organism] OR metagenomes[orgn]')
cont = Entrez.read(resp)
with open(snakemake.output[0], "w") as f:
    for item in cont["IdList"]:
        f.write("{}\n".format(item))