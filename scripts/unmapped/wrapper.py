
from snakemake.shell import shell

reformat_fasta_extra = snakemake.params.get("reformat_fasta_extra", "")

# Preprocessing command to run.
commands = [
            "reformat.sh in={snakemake.input} out={snakemake.output.fastq} unmappedonly primaryonly",
            "reformat.sh in={snakemake.output.fastq} out={snakemake.output.fasta} {reformat_fasta_extra}"
            ]

# Run preprocessing commands.
for cmd in commands:
  shell(cmd)
