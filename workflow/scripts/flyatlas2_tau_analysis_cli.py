#!/usr/bin/env python3.8

import click
import seaborn as sns
import matplotlib.pyplot as plt
from utils import flyatlas2_tau_analysis

"""
Downloads expression data from FlyAtlas2 'https://flyatlas.gla.ac.uk/FlyAtlas2/index.html?page=home'
"""


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--url',
              default="https://motif.mvls.gla.ac.uk/downloads/FlyAtlas2_gene_data.xlsx",
              help="URL link to FlyAtlas2 xlsx data")
@click.option('--sheet',
              default="FPKMs + SD",
              help="Sheet to keep from the Excel file")
def cli(url, sheet):
    "Command line interface"
    
    flyatlas2_tau_analysis(url, sheet)

if __name__ == "__main__":
    cli()
