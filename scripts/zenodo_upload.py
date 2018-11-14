
import os
import requests
import json

deposition_id = snakemake.params[0]

open_file = [open("{}".format(n), "rb") for n in snakemake.input]
files = list(zip(file, open_file))

r = requests.post('https://zenodo.org/api/deposit/depositions/{}/files'.format(deposition_id),
                  params = {'access_token': os.environ['ZENODO_PAT']},
                  files = files)
print(r.status_code)
