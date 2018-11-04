#!/usr/bin/env python3

import sys
from subprocess import Popen, PIPE
from common import filters

rm = "RepeatMasker -qq -pa {threads} {input} -dir {outdir}"
noreps = "if head -n 1 {output} | grep -q 'There were no repetitive sequences detected' then ln -sr {input} {masked} fi"

rm = rm.format(threads = str(snakemake.threads), input = snakemake.input["fa"], outdir = snakemake.params["outdir"])
noreps = noreps.format(output = snakemake.output["out"], input = snakemake.input["fa"], masked = snakemake.output["masked"])

p = Popen('/bin/bash', shell = False, universal_newlines = True, stdin = PIPE, stdout = PIPE, stderr = PIPE)
output, error = p.communicate(rm)
if p.returncode != 0:
   raise Exception("Job failed\n" + output + error)
else:
	print(output)

filters.run_filter_n(snakemake.input, snakemake.output, snakemake.params)
