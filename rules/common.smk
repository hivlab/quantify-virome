
# Load required libraries
import os.path
import json
from os import listdir
import glob
import pandas as pd
shell.executable("bash")

# Load configuration file with sample and path info
configfile: "config.yml"
samples = pd.read_table(config["samples"], sep = "\s+", index_col = "sample", dtype = str)
sample_ids = samples.index.values.tolist()
n_files = config["split_fasta"]["n_files"]

wildcard_constraints:
    n = "\d+"

# Create logs dir
import os
if not os.path.exists("logs"):
    os.makedirs("logs")
