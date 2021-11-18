#!/usr/bin/env python3

"Filter the orthogroups containing cell adhesion molecules (CAMs), as described in FlyXCDB"

import pandas as pd

def get_cams_orthgroups(cams, orthogroups):
    ""

    # Read the surface proteins downloaded from FlyXCDB
    cams = pd.read_csv(cams)

    # Filter the first 320, which are likely CAMs
    cams = cams.iloc[:320, 2]

    # Load the orthogroups obtained from Orthofinder
    ort_grps = pd.read_csv(orthogroups, sep='\t')

    # Keep only the orthogroups were Drosophila melanogaster is present
    ort_grps = ort_grps.loc[:, ['Orthogroup', 'Drosophila_melanogaster.BDGP6.32.pep.all']].dropna()

    # Store the orthogroup names into a Series
    dmel_ort_grps = ort_grps['Orthogroup']

    # Store FlyBase identifiers for each orthogroup
    flybase_ids = [id.split('|')[3].lstrip('gene_')
                   for id in ort_grps['Drosophila_melanogaster.BDGP6.32.pep.all']]

    # Create a dictionary mapping orthogroup name and D. melanogaster gene symbol
    ort_dict = {key: value
                for key, value in zip(flybase_ids, dmel_ort_grps)}
    # Open a new file to store CAMs orthogroups
    with open("cams_orthogroups.csv", 'w') as fh:
        # For each CAM, write the orthogroup name
        for cam in cams:
            if cam in ort_dict.keys():
                fh.write(f'{ort_dict[cam]}\n')

get_cams_orthgroups('~/flyxcdb.csv', 'Orthogroups/Orthogroups.tsv')
