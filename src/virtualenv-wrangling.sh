#!/bin/bash

module load python-3.6.0
virtualenv env
source env/bin/activate
pip install snakemake
deactivate
module purge
