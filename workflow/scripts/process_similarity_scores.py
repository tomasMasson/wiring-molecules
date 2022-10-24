#!/usr/bin/env python3.8

"""
Process the similarity scores from a Diamond blastp and returns a matrix with
standarized scores
"""

import click
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def process_similarity_scores(matrix, output):
    """
    Transforms a tab-separated Diamond blastp results file into a matrix with
    normalized values.
    """

    # Define columns names
    columns = ["qseqid", "sseqid", "pident", "length",
               "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore"]
    # Read data into a dataframe
    df = pd.read_csv(matrix,
                     sep="\t",
                     names=columns,
                     usecols=[0, 1, 11])
    # Extract the species name from query identifiers
    df["Species"] = df["qseqid"].str.split("_", expand=True).loc[:, 1]
    # Pivot dataframe to separate different specie scores
    df = df.pivot_table(index='sseqid', columns="Species", values="bitscore")
    # Set query protein name as index
    indexes = [index[2] for index in df.index.str.split("|")]
    df.index = indexes
    # Fill NaN values
    df = df.fillna(20.4)
    # Normalize blastp scores to Drosophila melanogaster values
    gene_norm = df["DROME"]
    df_n = df.divide(gene_norm, axis="index")
    # Normalize for species divergence (Z-scaling across columns)
    df_n = df_n.subtract(df_n.mean(axis=0), axis=1).divide(df_n.std(axis=0), axis=1)
    # Save matrix to output file
    df_n = df_n.drop('DROME', axis=1)
    df_n = df_n.transpose()
    df_n.to_csv(output)
    # print("Plotting correlations")
    # corr = df_n.corr(method='pearson')
    # corr.to_csv(output)
    # p = sns.clustermap(df_n)
    # fig = p.get_figure()
    # p.savefig("tmp.png")


@click.command()
@click.option("-i",
              "--input",
              help="Input file with Diamond blastp results.")
@click.option("-o",
              "--output",
              help="Output file")
def cli(input, output):
    "Command line interface"
    process_similarity_scores(input, output)


if __name__ == "__main__":
    cli()
