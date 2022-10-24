#!/usr/bin/env python3

import click
import pandas as pd
from functools import reduce


def build_conservation_matrix(input, output):
    """
    Starts from an aggregated Diamond tsv file
    and generates a pandas-like matrix where
    each row is reference proteome protein and
    columns store different species

    Only E-values lower than 0.001 are deemed
    significant
    """

    # Define columns names
    columns = ["qseqid", "sseqid", "pident", "length",
               "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore"]
    # Read data into a dataframe
    df = pd.read_csv(input,
                     sep="\t",
                     names=columns,
                     usecols=[0, 1, 10, 11])
    # Drop non-significant E-values
    df.loc[df.evalue > 0.00001, "bitscore"] = 0
    # Store species names into a column and a set
    df["species"] = df["sseqid"].str.split("_", expand=True)[1]
    species = set(df["sseqid"].str.split("_", expand=True)[1])
    # Split data coming from different proteomes
    data_col = []
    for s in species:
        dff = df[df["species"] == s]
        dff.rename(columns={"bitscore": s}, inplace=True)
        data_col.append(dff.loc[:, ["qseqid", s]])
    # Aggregate data into a single matrix
    # Each species is a single column
    matrix = reduce(lambda x, y: pd.merge(x, y, on="qseqid", how="outer"), data_col).fillna(0)
    # Save matrix to output file
    matrix.to_csv(output, index=False)


@click.command()
@click.option("-i",
              "--input",
              help="Input file with Diamond blastp results.")
@click.option("-o",
              "--output",
              help="Output file name.")
def cli(input, output):
    "Command line interface"
    build_conservation_matrix(input, output)


if __name__ == "__main__":
    cli()
