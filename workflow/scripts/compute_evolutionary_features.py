#!/usr/bin/env python3

import click
import pandas as pd
from collections import Counter


def normalize(series):
    """
    Min-Max normalization
    """

    max_value = series.max()
    min_value = series.min()
    norm_series = (series - min_value) / (max_value - min_value)

    return norm_series


def compute_size(row):
    """
    Returns the size of a protein list stored on a DataFrame column
    """

    # Extract ortholog proteins columns
    n_proteins = len(
                     row["List_orthologs"].split(",")
                   )

    return n_proteins


def compute_retention(row):

    proteins = row["List_orthologs"].split(",")
    species = len(set(
              [protein.split("_0:")[0]
               for protein in proteins]
              ))
    return species


def compute_duplicability(row):

    proteins = row["List_orthologs"].split(",")
    species = Counter(
              [protein.split("_0:")[0]
               for protein in proteins]
              )
    dups_species = len([key
                        for key in species
                        if species[key] > 1])

    return dups_species


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input',
              help="Orthogroups input data")
# CLI main function
def cli(input):
    pass

    # Load data
    df = pd.read_csv(input, sep="\t")
    # Calculate evolutionary features
    df["Orthogroup_size"] = normalize(df.apply(compute_size, axis=1))
    df["Retention"] = normalize(df.apply(compute_retention, axis=1))
    df["Duplicability"] = normalize(df.apply(compute_duplicability, axis=1))

    df.to_csv("evolutionary_features.csv", index=False)


if __name__ == "__main__":
    cli()
