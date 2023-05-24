#!/usr/bin/env python3

"Download genome datasets from from NCBI FTP"

import requests
import click
import pandas as pd


def download_url(url, output):
    """
    Downloads URL content into an output file
    """

    # Get a response from the URL
    response = requests.get(url)

    # Write the response content into the output file
    with open(output, 'wb') as fh:
        fh.write(response.content)


def download_ncbi_genomes(url_list):
    """
    Downloads NCBI entries based on a URLs list
    """

    # Load URLs list
    df = pd.read_csv(url_list, names=["organism", "url"])
    # Download single endpoint
    for row in df.iterrows():
        url = row[1].url
        output = f"{row[1].organism}.fa.gz"
        print(f"Downloading '{output}'\n---")
        download_url(url, output)
        print(f"Finished downloading '{output}'\n---")


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-l", "--list",
              help="URLs list")
# CLI main function
def cli(list):
    """
    Downloads genome dataset from NCBI FTP based on a list of URLs

    Invocation example:

    ./download_ncbi_genomes.py --list diptera_urls.csv --folder ../results/
    """

    download_ncbi_genomes(list)


if __name__ == '__main__':
    cli()
