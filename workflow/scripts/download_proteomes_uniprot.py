#!/usr/bin/env python3

"Download proteome data from UniProtKB"

import urllib.request
import click


# Set CLI parameters
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-s", "--species", help="Species list to retrieve from UniProtKB")
def download_uniprot_proteomes(species):
    """
    Download a list of input proteomes from UNIPROT ftp server
    """

    # Define the path for Uniprot database
    ftp = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"

    # Parse proteomes list
    with open(species, 'r') as fh:
        # Store data into a list
        proteomes = [(line.split("\t")[0], line.split("\t")[2].rstrip())
                     for line in fh]
    # Iterate over data rows
    for item in proteomes[1:]:
        # Compose FTP endpoint name
        endpoint = item[0] + "_" + item[1] + ".fasta.gz"
        # Compose URL
        url = ftp + item[0] + "/" + endpoint
        try:
            # Make request
            response = urllib.request.urlopen(url)
            # Read request content
            zipcontent = response.read()
            print(f"Currently downloading {endpoint}")
            # Write request content to file
            with open(endpoint, 'wb') as fd:
                fd.write(zipcontent)
            print(f"Just finished with {endpoint}")
        except urllib.error.HTTPError as e:
            print(e.reason)


# Execute CLI command
if __name__ == '__main__':
    download_uniprot_proteomes()
