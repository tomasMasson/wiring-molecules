#!/usr/bin/env python3


import click
import pandas as pd


def get_consensus_surfaceomes(data):
    """
    Reads the output table from ratiometric_classifier.py and returns a table
    with the consensus (present in all four controls) proteome for each sample
    """

    # Ingest protein table
    df = pd.read_csv(data,
                     names=["Protein", "Sample",
                            "Control"])
    # Group proteins by sample
    grps = df.groupby("Sample")
    # Initialize the dict to store individual proteomes as lists
    protein_lists = {}
    # Iterate over sample lists to assess protein repoducibility among controls
    for grp in grps:
        # Store sample name
        sample = grp[0]
        # Store protein counts among control conditions
        counts = grp[1]["Protein"].value_counts()
        # Keep only proteins found in all four controls
        protein_lists[sample] = list(counts[counts == 4].index)
        # For sample T4_2, be least stringent and keep all proteins identified
        if sample == "t4_2":
            protein_lists[sample] = list(counts[counts >= 1].index)

    # Intersect proteins list from replicates to obtain the consensus table
    t4_prot = set(protein_lists["t4_1"]).intersection(protein_lists["t4_2"])
    t5_prot = set(protein_lists["t5_1"]).intersection(protein_lists["t5_2"])
    # Save final tables to output files (one for each type, T4 and T5)
    (pd.Series(list(t4_prot), name="Protein")
     .to_csv("t4_consensus_surfaceome.csv", index=False)
     )
    (pd.Series(list(t5_prot), name="Protein")
     .to_csv("t5_consensus_surfaceome.csv", index=False)
     )


@click.command()
@click.option("-i",
              "--input",
              help="Proteins table obtained after ratiometic classification")
def cli(input):
    "Command line interface"
    get_consensus_surfaceomes(input)


if __name__ == "__main__":
    cli()
