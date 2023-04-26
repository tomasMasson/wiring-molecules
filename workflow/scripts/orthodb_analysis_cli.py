#!/usr/bin/env python3.8

"Download, process and analyze orthogroups data from OrthoDB"

import click
import pandas as pd
from orthodb_utils import fetch_orthodb_data


# Command line interface
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-p",
              "--pipeline",
              type=str,
              help="Stages of the pipeline to be executed [full, download, analysis]")
@click.option("-u",
              "--url",
              type=str,
              default="https://data.orthodb.org/download/",
              help="OrthoDB URL that will be used to download orthogroups data")
@click.option("-o",
              "--outdir",
              type=str,
              default="../resources/orthodb",
              help="Directory address to save OrthoDB data")
def cli(pipeline: str, url: str, outdir:str) -> None:
    "OrthoDB analysis pipeline"
    
    if not pipeline:
        print("Please specify the pipeline step that you would like to run (full, download or analysis)")
        return None
    elif pipeline == "full":
        print("Full pipeline will be executed")
        fetch_orthodb_data(url, outdir)
        print("Data analysis")
    elif pipeline == "download":
        print("OrthoDB data will be downloaded")
        fetch_orthodb_data(url, outdir)
    elif pipeline == "analysis":
        print("Data analysis")
    return None


if __name__ == "__main__":
    cli()
