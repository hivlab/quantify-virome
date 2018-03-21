#!/usr/bin/env python3
"""
Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --cluster slurm_scheduler.py
"""

import os
import sys
import warnings


def eprint(*args, **kwargs):
    print(*args, file = sys.stderr, **kwargs)

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]

job_properties = read_job_properties(jobscript)

cluster_param = {}

job_resources = job_properties["resources"]


if not "mem" in job_resources:
    warnings.warn("Rule {rule} has no memory specified, set to default.".format(**job_properties))

# do something useful with the threads
cluster_param["threads"] = job_properties.get("threads",1)
cluster_param['days'] = job_resources.get("days", job_properties["cluster"]['days'])
cluster_param['hours'] = job_resources.get("hours", job_properties["cluster"]['hours'])
cluster_param['mem'] = int(job_resources.get("mem", 10)) + 5 #GB + overhead
cluster_param['name'] = job_properties['rule']

cluster_param['profile'] = job_properties["cluster"]['profile']

# access property defined in the cluster configuration file (Snakemake >=3.6.0)
#job_properties["cluster"]["profile"]

eprint("Submit job with parameters:\n"+"\n".join(["\t{} : {}".format(key,cluster_param[key]) for key in cluster_param]))

os.system("sbatch -p {profile} --parsable -c {threads} --time={days:d}-{hours:02d}:00:00 --mem={mem}g --job-name={name} {script}".format(script=jobscript, **cluster_param))
#
