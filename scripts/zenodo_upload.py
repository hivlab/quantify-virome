
import os
import requests
import json

deposition_id = snakemake.params[0]

open_file = [open("{}".format(n), "rb") for n in snakemake.input]
files = list(zip(snakemake.input, open_file))

r = requests.get('https://zenodo.org/api/deposit/depositions', params={'access_token': os.environ['ZENODO_PAT']})
r.status_code
r.json()

headers = {"Content-Type": "application/json"}
r = requests.post('https://zenodo.org/api/deposit/depositions', params={'access_token': os.environ['ZENODO_PAT']}, json={}, headers=headers)
r.status_code
r.json()

data = {}
files = {'file': open('results/SRR5580357_phages_1.csv', 'rb')}
r = requests.post('https://zenodo.org/api/deposit/depositions/{}/files'.format(deposition_id),
                  params = {'access_token': os.environ['ZENODO_PAT']},
                  data = data,
                  files = files)
print(r.status_code)
r = requests.post('https://zenodo.org/api/deposit/depositions/{}/actions/publish'.format(deposition_id),
                      params={'access_token': os.environ['ZENODO_PAT']})
print(r.status_code)
