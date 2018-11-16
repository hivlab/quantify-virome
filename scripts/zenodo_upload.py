import os
import requests
import gzip
import re
from subprocess import Popen, PIPE
import hashlib

if not os.environ['ZENODO_PAT']:
      raise ValueError("Missing ZENODO_PAT environment variable with zenodo API access token!")

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

# Upload depository id
deposition_id = snakemake.params[0]

# Create tar.gz file for upload
files = snakemake.input
zipfile = list(set([re.sub("_\d+", "", file) for file in files]))[0] + ".tar.gz"
cmd = ['tar', '-cvzf', zipfile] + files
p = Popen(cmd, stdout = PIPE, stderr = PIPE)
stout, stderr = p.communicate()

# Check if file is present
# calculate md5 checksum for local file
hash = md5(zipfile)

# Compose files query and upload url
base_url = 'https://zenodo.org/api'
url = os.path.join(base_url, 'deposit/depositions/{}/files'.format(deposition_id))

# Setup access token
params = {'access_token': os.environ['ZENODO_PAT']}

# Get info for remote files
r = requests.get(url, params = params)
filename = [deposit['filename'] for deposit in r.json()]

# Upload, if file is not present
if os.path.basename(zipfile) not in filename:
    with open(zipfile, "rb") as handle:
        r = requests.post(url, params = params,
                          data = {'filename': str(zipfile)},
                          files = {'file': handle})

        if r.status_code != 201:
                raise requests.HTTPError(f"Error in data upload, status code: {r.status_code}   {r.json()['message']}")

else:
    print("Doing nothing. File {} is already uploaded!\nPlease delete local and remote copy of the file\nif you wish to upload new version.".format(os.path.basename(zipfile)))
