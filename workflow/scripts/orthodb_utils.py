import shutil
import os
from pathlib import Path
from pprint import pprint

import requests
import pandas as pd
from bs4 import BeautifulSoup


def get_url_links(url: str) -> list:
    "Returns all the URL links (markerd with an <a> tag) inside HTMl content"
    # Load the HTML file
    with requests.get(url) as r:
        soup = BeautifulSoup(r.content, "html.parser")
    # Extract all hyperlinks
    links = [link.contents[0]
             for link in soup.find_all("a")
             # Ignore large fasta files from sequences
             if not link.contents[0].endswith("fasta.tab.gz")]
    return links


def check_orthodb_exists(outdir:str) -> bool:
    "Assess if OrthoDB data already exists or should downloaded"
    # Check if the output directory exist in Resources
    p = Path(outdir)
    if p.exists():
        print(f"{outdir} already exists")
        print("Data will not be downloaded again")
        return True
    else:
        print("Orthodb data download will start")
        print(f"{outdir} has been created")
        os.mkdir(outdir)
        return False


def download_compressed_file(url: str, outfile: str) -> None:
    "Download OrthoDB data"
    # Connect to url data stream
    with requests.get(url, stream=True) as r:
        # Save data to output file as chunks
        with open(outfile, "wb") as fh:
            for chunk in r.iter_content(chunk_size=8192):
                fh.write(chunk)
    return None


def move_files(files: list, outdir: str) -> None:
    "Move downloaded files into the output directory"
    p = Path(outdir)
    for file in files:
        outfile = p / Path(file)
        shutil.move(file, outdir)
        print(f"{outfile} was moved into {outdir}")
    return None


def fetch_orthodb_data(url: str, outdir: str) -> None:
    "Download and stored files from OrthoDB server into output folder"
    # Check if an updated URL has been provided
    if not url:
        url = "https://data.orthodb.org/download/"
    # Check if files have been already downloaded
    if not check_orthodb_exists(outdir):
        # Get the names for the link
        odb_links = get_url_links(url)
        # Download actual files
        print("Files to be downloaded:")
        pprint(odb_links)
        for link in odb_links:
            print(f"Downloading {link}")
            endpoint = f"{url}{link}"
            download_compressed_file(endpoint, link)
        # Move everything into the desired folder
        move_files(odb_links, outdir)
    return None


def process_orthodb_data(data):
    pass


def filter_orthogroups_table(data, taxid):
    """
    Returns a filtered DataFrame from OrthoDB orthogroups table (odb10v1_OGs)
    based on the TaxID provided
    """

    # Column names
    COLUMNS = ["OG", "Level", "OG_name"]
    # Load orthogroups table
    df = pd.read_csv(data,
                     sep="\t",
                     names=COLUMNS)
    # Filter orthogroups based on Taxid
    fdf = df[df["Level"] == taxid]
    print(fdf)
    return fdf
