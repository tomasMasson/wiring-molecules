#!/usr/bin/env python3

"Download proteome data from Ensembl Genomes"

import urllib.request
import requests
import click

# Set CLI parameters
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-s", "--species", help="Species list to retrieve from UniProtKB")

def download_ensembl_proteomes(species):
    """
    Download a list of input proteomes from Ensembl ftp server
    """

    # Define the path for Uniprot database
    ftp = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/metazoa/fasta/"

    # Parse proteomes list
    with open(species, 'r') as fh:
        # Store data into a list
        proteomes = [line.split('\t')[1]
                     for line in fh]
    # Iterate over data rows
    for item in proteomes[1:]:
        # Compose FTP endpoint name
        folder = item + '/pep/'
        # Compose folder URL
        url = ftp + folder
        # Make folder request
        response = requests.get(url)
        # Read folder request content
        for line in response.text.splitlines():
            if '.fa.gz' in line:
                endpoint = line.split('"')[1]
        # Compose endpoint URL
        url = ftp + folder + endpoint
        # Make folder request
        response = urllib.request.urlopen(url)
        print(f"Currently downloading {endpoint}")
        # Read endpoint request content
        zipcontent = response.read()
        # Write request content to file
        with open(endpoint, 'wb') as fd:
            fd.write(zipcontent)
        print(f"Just finished with {endpoint}")

# Execute CLI command
if __name__ == '__main__':
    download_ensembl_proteomes()
