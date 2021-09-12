#!/usr/bin/env python3

"Module docstring"

from pathlib import Path
import gzip

p = Path("data/uniprot_proteomes/")
for child in p.iterdir():
    proteome = (str(child).split("/")[2].split('.gz')[0])
    with gzip.open(child, 'rb') as f_in:
        content = f_in.read()
        with open(proteome, 'wb') as f_out:
            f_out.write(content)
