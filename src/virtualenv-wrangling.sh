#!/bin/bash

module load python-3.6.0
virtualenv env
env/bin/pip install snakemake
env/bin/pip freeze > requirements.txt
module purge

