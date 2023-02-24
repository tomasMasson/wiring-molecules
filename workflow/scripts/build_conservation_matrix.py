#!/usr/bin/env python3

import click
import pandas as pd
import numpy as np
from functools import reduce


def get_mapping(mapping):
    """
    Load mapping between UniProt and FlyBase accessions
    """

    df = pd.read_csv(mapping,
                     skiprows=19, sep="\t",
                     names=["Symbol", "Organism",
                            "FBgn", "UniProt"]
                     )
    dff = df.dropna().iloc[:, 2:]
    return dff


def build_conservation_matrix(input, mapping, reference, output):
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
    df["species"] = df["sseqid"].str.split("|", expand=True)[3]
    species = set(df["sseqid"].str.split("|", expand=True)[3])
    # Split data coming from different proteomes
    data_col = []
    for s in species:
        dff = df[df["species"] == s]
        dff.rename(columns={"bitscore": s}, inplace=True)
        data_col.append(dff.loc[:, ["qseqid", s]])
    maps = get_mapping(mapping)
    # Aggregate data into a single matrix
    # Each species is a single column
    matrix = reduce(lambda x, y: pd.merge(x, y, on="qseqid", how="outer"), data_col).fillna(0)
    # Define protein names column
    matrix.insert(0, "Protein_ID", matrix.qseqid.str.split("|", expand=True)[1])
    matrix = matrix.drop("qseqid", axis=1)
    # Map UniProt accessions into FlyBase FBgn identifiers
    genes = dict(zip(maps["UniProt"], maps["FBgn"]))
    matrix.replace(genes)
    # Set protein names as indexes
    matrix.set_index("Protein_ID", inplace=True)
    # Fill missing values
    matrix.fillna(1, inplace=True)
    # Compute relative bitscores againts a reference organism
    # (each protein/row has a different length, biasing the score comparison)
    gene_norm = matrix[reference]
    nmatrix = np.log2(matrix.divide(gene_norm, axis="index"))
    # Normalize each species/column to account for their different phylogenetic histories)
    nmatrix = (nmatrix - nmatrix.mean()) / nmatrix.std()
    # Drop reference column used for normalization
    nmatrix.dropn(columns=[reference], inplace=True)
    # Save matrix and correlations into output file
    nmatrix.to_csv(output)
    corr = nmatrix.T.corr()
    corr.to_csv("correlations.csv")


@click.command()
@click.option("-i",
              "--input",
              help="Input file with Diamond blastp results.")
@click.option("-m",
              "--mapping",
              help="Identifiers mapping: UniProt -> FlyBase.")
@click.option("-r",
              "--reference",
              help="Reference species for normalization.")
@click.option("-o",
              "--output",
              help="Output file name.")
def cli(input, mapping, reference, output):
    "Command line interface"
    build_conservation_matrix(input, mapping, reference, output)


if __name__ == "__main__":
    cli()
