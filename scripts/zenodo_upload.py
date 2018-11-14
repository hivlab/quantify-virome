
import os
import requests
import json

headers = {"Content-Type": "application/json"}
r = requests.get('https://zenodo.org/api/deposit/depositions',
                  params = {'access_token': os.environ['ZENODO_PAT']},
                  json = {}, headers = headers)

deposition_id = snakemake.params[0]
r = requests.post('https://zenodo.org/api/deposit/depositions/{}/files'.format(deposition_id),
                  params = {'access_token': os.environ['ZENODO_PAT']},
                  files = files)
