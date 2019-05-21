
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout = False, stderr = True)

shell("(cd-hit-est"
      " -i {snakemake.input}"
      " -o {snakemake.output.repres}"
      " -T {snakemake.threads}"
      " {extra}) {log}")
