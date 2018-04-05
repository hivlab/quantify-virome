#!/bin/bash
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake
module load python-3.6.0
virtualenv env
source env/bin/activate
pip install cookiecutter
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
deactivate
