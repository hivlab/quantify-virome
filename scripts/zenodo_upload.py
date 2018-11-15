import os
import requests
import gzip
import re
from subprocess import Popen
import hashlib

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
p = Popen(cmd)

# Check if file is present
# calculate md5 checksum for local file
hash = md5(zipfile)

# Get info for remote files
r = requests.get('https://zenodo.org/api/deposit/depositions/{}/files'.format(deposition_id),
                        params={'access_token': os.environ['ZENODO_PAT']})
fchk = [{deposit['filename'], deposit['checksum']} for deposit in r.json()]
filename, checksum = zip(*fchk)

if hash not in list(checksum):
    with open(zipfile, "rb") as handle:
        r = requests.post('https://zenodo.org/api/deposit/depositions/{}/files'.format(deposition_id),
                          params = {'access_token': os.environ['ZENODO_PAT']},
                          data = {'filename': str(zipfile)},
                          files = {'file': handle})

        if r.status_code != 201:
                raise requests.HTTPError(f"Error in data upload, status code: {r.status_code}   {r.json()['message']}")

else:
    print("File {} is already uploaded as {}".format(zipfile, [f for f in list(filename) if f == os.path.basename(list(filename))][0]))
