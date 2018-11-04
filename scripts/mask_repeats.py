#!/usr/bin/env python3

from subprocess import Popen, PIPE
from common import filters

rm = "RepeatMasker -qq -pa {threads} {input} -dir {outdir}"
rm = rm.format(threads = str(snakemake.threads), input = snakemake.input["fa"], outdir = snakemake.params["outdir"])

# Run repeatmasker
p = Popen(rm.split(" "), stdout = PIPE, stderr = PIPE, universal_newlines=True)
output, error = p.communicate()
if p.returncode != 0:
   raise Exception("Job failed\n" + output + error)
else:
  print(output)

# Filter repeatmasker output
print("Filtering runs\n")
filters.run_filter_n(snakemake.input, snakemake.output, snakemake.params)
