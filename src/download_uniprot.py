#!/usr/bin/env python3

"Module docstring"

import urllib.request

# Define the path for Uniprot database

UNIPROT_FTP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"

# Import the list of proteomes
with open("data/reference_proteomes_arthopoda_uniprot.tab", "r") as fh:
    proteomes = [(line.split("\t")[0], line.split("\t")[2])
                 for line in fh]

for item in proteomes[1:]:
    endpoint = item[0] + "_" + item[1] + ".fasta.gz"
    url = UNIPROT_FTP + item[0] + "/" + endpoint
    response = urllib.request.urlopen(url)
    zipcontent= response.read()
    print(f"Currently downloading {endpoint}")
    with open(endpoint, 'wb') as fd:
        fd.write(zipcontent)
    print(f"Just finished with {endpoint}")
