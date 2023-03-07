#!/usr/bin/env python3

"Download full proteome data from Ensembl Genomes"

import urllib.request
import requests
import click
import pandas as pd
from bs4 import BeautifulSoup

def select_ensembl_species(species, table):
    """
    Filters the species of interests, called group, from a table containing all the species in the Current release
    """

    # Read the species table from Ensembl Genomes
    df = pd.read_csv(table, sep='\t', index_col=False)

    # Filter out the species that are not in the group list
    species_filt = df[df.species.isin(species)].species

    return species_filt


def download_file(url, output):
    """
    Downloads the content from a URL into an output file
    """

    # Get a response from the URL
    response = requests.get(url)

    # Write the response content into the output file
    with open(output, 'w') as fh:
        fh.write(response.text)


def download_ensembl_proteome(organism):
    """
    Downloads a single Ensembl Genomes entry corresponding to the organism input. File is save into output.
    """

    # Ensembl Genomes fasta folder
    ftp = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/current/fasta/"

    # Compose FTP folder name for the organism
    folder_url = ftp + organism + '/pep/'

    # Make folder request
    response = requests.get(folder_url)

    # Read folder request content
    for line in response.text.splitlines():
        # Check for the compressed sequence file
        if '.fa.gz' in line:
            endpoint = line.split('"')[1]

    # Compose endpoint URL
    url = folder_url + endpoint

    # Make folder request
    response = urllib.request.urlopen(url)
    print(f"Downloading '{endpoint}'\n---")

    # Read endpoint request content
    zipcontent = response.read()

    # Write request content to file
    output = f"{organism}.fa.gz"
    with open(output, 'wb') as fd:
        fd.write(zipcontent)
    print(f"Finished downloading {output}\n---")


def download_ensembl_species_list(output):
    """
    Download the species list available in the Current release at Ensembl ftp server
    """

    # Path to the current release of Ensembl Genomes
    ftp = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/current/species_EnsemblMetazoa.txt"

    # Download the species list from the Current release
    download_file(ftp, output)

    # Print where the content was save
    print(f"Ensembl Genomes table was saved to '{output}'\n---")


def download_ensembl_current_proteomes(species, ensembl_list):
    """
    Download a list of input proteomes from the Ensembl ftp server
    """

    # Read the species input into a list
    with open(species, 'r') as fh:
        species_list = [line.split()[0] for line in fh]

    # Filter species of interest
    species_paths = select_ensembl_species(species_list, ensembl_list)

    # Iterate over species data
    for path in species_paths:
        download_ensembl_proteome(path)

    print('Data download is finished')


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-e", "--ensembl_list",
              help="Name for the species list retrieved from Ensembl Genomes")
@click.option("-s", "--species",
              help="Species list used to download proteome data from Ensembl Genomes")

# CLI main function
def command_line_interface(ensembl_list, species):
    """
    Downloads whole-proteome data from Ensembl Genomes based on a list of input identifiers.

    Invocation example:

    ./download_ensembl_proteomes.py --ensembl-list ensembl_species_list.tsv --species species_list.csv
    """

    download_ensembl_species_list(ensembl_list)

    download_ensembl_current_proteomes(species, ensembl_list)

if __name__ == '__main__':
    command_line_interface()
