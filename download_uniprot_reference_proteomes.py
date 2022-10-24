#!/usr/bin/env python3

import requests
from bs4 import BeautifulSoup

url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
r = requests.get(url)
content = r.content
soup = BeautifulSoup(content, "html.parser")
results = soup.find_all("a")
for i in results[5:]:
    print(str(i).split('"')[1][:-1])
