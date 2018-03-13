#!/bin/bash

module load python-3.6.0
virtualenv -p python3 ~/Projects/vs/.venv
source ~/Projects/vs/.venv/bin/activate
pip install snakemake
pip install fontconfig


##to deactivate and remove
cd ~/Projects/vs/
deactivate
rm -r ~/Projects/vs/.venv
module purge
