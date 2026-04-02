#!/usr/bin/env python3

"""
Extract the members of each orthogroup stored at OrthoDB for a given taxon
"""

import click
import subprocess
from numpy import add
import pandas as pd


def add_orthodb_species(df: pd.DataFrame) -> pd.DataFrame:
    "Adds OrthoDB species id as column, as parsed from gene identifier"
    return (df
            .assign(species=df.gene.str.split(":", expand=True)[0])
           )


def group_orthogroups(df: pd.DataFrame, column: str) -> pd.DataFrame:
    # Group genes by orthogroup
    groups = df.groupby(column)
    data = []
    COLUMNS = ["og", "species"]
    for group in groups:
        # Store species name
        og = group[0]
        # Store gene name
        species = ",".join(list(group[1]["species"]))
        # Store number of genes in the orthogroup
        data.append([og, species])
    dff = pd.DataFrame(data, columns=COLUMNS)
    return dff


def group_orthogroups2(df: pd.DataFrame, column: str) -> pd.DataFrame:
    # Group genes by orthogroup
    groups = df.groupby(column)
    data = []
    COLUMNS = ["og", "species"]
    for group in groups:
        # Store species name
        og = group[0]
        # Store gene name
        species = ",".join(list(group[1]["og"]))
        # Store number of genes in the orthogroup
        data.append([og, species])
    dff = pd.DataFrame(data, columns=COLUMNS)
    return dff


def map_fbgn(df, idx, uniprot, fbgn):
    idx = pd.read_csv(idx, sep="\t", names=["OG", "Gene"])
    uniprot = pd.read_csv(uniprot, sep="\t", names=["Gene", "UniProt", "DB"])
    fbgn = pd.read_csv(fbgn, sep="\t", skiprows=5, usecols=[2, 5], names=["FBgn", "UniProt"]).dropna()
    mappings = pd.merge(pd.merge(idx, uniprot), fbgn)
    return df.replace(dict(zip(mappings.OG, mappings.FBgn)))


def preprocess_dataframe(data: str, labels: str, uniprot:str, fbgn:str, outfile: str) -> None:

    df = pd.read_csv(data,
                     sep="\t",
                     names=["og", "gene"])
    return (df
            # .pipe(add_orthodb_species)
            # .pipe(group_orthogroups, "og")
            .pipe(map_fbgn, labels, uniprot, fbgn)
            .to_csv(outfile, index=False)
            )


# def preprocess_dataframe2(data: str, labels: str, uniprot:str, fbgn:str, outfile: str) -> None:

#     df = pd.read_csv(data,
#                      sep="\t",
#                      names=["og", "gene"])
#     return (df
#             .pipe(add_orthodb_species)
#             .pipe(group_orthogroups2, "species")
#             .pipe(map_fbgn, labels, uniprot, fbgn)
#             .to_csv(f"{outfile}2", index=False)
#             )


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-og2g', '--og2genes',
              help="Mapping from orthogroups to gene names")
@click.option('-l', '--labels')
@click.option('-u', '--uniprot')
@click.option('-f', '--fbgn')
@click.option('-o', '--outfile',
              help="Output file name")
def cli(og2genes, labels, uniprot, fbgn, outfile):
    """
    Computes proteins members of each orthogroup at OrthoDB for a given taxon
    """

    preprocess_dataframe(og2genes, labels, uniprot, fbgn, outfile)
    # preprocess_dataframe2(og2genes, labels, uniprot, fbgn, outfile)


if __name__ == '__main__':
    cli()
