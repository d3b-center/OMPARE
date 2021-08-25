# This script uses the id from the maf_manifest.tsv and query the GDC API data endpoint for download
# the downlaoded the files will be in their respective folder - needed to move them out of their own
# folder and put in `maf_files`

import requests
import json
import pandas as pd
import re
import os

# read in the manifest as a dataframe and subset the id column
maf_manifest = pd.read_csv('../results/tcga_not_in_pbta_maf_manifest.tsv', sep='\t')
ids = maf_manifest["id"].values.tolist()
# find the unique ids for the next step
ids = list(set(ids))

data_endpt = "https://api.gdc.cancer.gov/data"
params = {"ids": ids}
response = requests.post(data_endpt,
                        data = json.dumps(params),
                        headers={
                            "Content-Type": "application/json"
                            })

response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

path_to_file = os.path.join("../../scratch", file_name)
with open(path_to_file, "wb") as output_file:
    output_file.write(response.content)
