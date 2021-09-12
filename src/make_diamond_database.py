#!/usr/bin/env python3

"Module docstring"

from pathlib import Path
import subprocess

p = Path("phylogenetic_profiles/proteomes/")
for child in p.iterdir():
    database = "phylogenetic_profiles/databases/" + str(child).split("/")[2]
#    # Build Diamond databases
#    build_db = f"diamond makedb --in {child} --db {database}"
#    subprocess.run(build_db, shell=True, check=True)
    score = "phylogenetic_profiles/scores/" + str(child).split("/")[2] + ".score"
    search = f"diamond blastp --query phylogenetic_profiles/UP000000803_7227.fasta --db {database} --max-target-seqs 1 --outfmt 6 --out {score}"
    subprocess.run(search, shell=True, check=True)
