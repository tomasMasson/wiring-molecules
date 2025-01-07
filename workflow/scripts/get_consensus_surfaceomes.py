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
                     names=["Protein",
                            "signal_ratio",
                            "Sample",
                            "Control"])
    t4 = df[df.Sample.isin(["t4_1", "t4_2"])]
    m = [i[1].signal_ratio.median() for i in t4.groupby("Protein")]
    p = [i[0] for i in t4.groupby("Protein")]
    t4_dic = dict(zip(p, m))
    t5 = df[df.Sample.isin(["t5_1", "t5_2"])]
    m = [i[1].signal_ratio.median() for i in t5.groupby("Protein")]
    p = [i[0] for i in t5.groupby("Protein")]
    t5_dic = dict(zip(p, m))
    # Group proteins by sample
    grps = df.groupby("Sample")
    # Initialize the dict to store individual proteomes as lists
    protein_lists = {}
    # Iterate over sample lists to assess protein reproducibility among controls
    for grp in grps:
        # Store sample name
        sample = grp[0]
        # Store protein counts among control conditions
        counts = grp[1]["Protein"].value_counts()
        # Keep only proteins found in all four controls
        protein_lists[sample] = list(counts[counts == 4].index)
        # # For sample T4_2, be least stringent and keep all proteins identified
        if sample == "t4_2":
            protein_lists[sample] = list(counts[counts >= 1].index)

    # Intersect proteins list from replicates to obtain the consensus table
    t4_prot = set(protein_lists["t4_1"]).intersection(protein_lists["t4_2"])
    t5_prot = set(protein_lists["t5_1"]).intersection(protein_lists["t5_2"])
    # # Save final tables to output files (one for each type, T4 and T5)

    t4 = pd.DataFrame(data=[t4_prot, t4_prot]).T
    t4.columns = ["Protein", "T4_signal"]
    t4.T4_signal.replace(t4_dic, inplace=True)
    t4 = t4.sort_values(by=["T4_signal"], ascending=False)
    t4.to_csv("t4_consensus_surfaceome.csv")

    t5 = pd.DataFrame(data=[t5_prot, t5_prot]).T
    t5.columns = ["Protein", "T5_signal"]
    t5.T5_signal.replace(t5_dic, inplace=True)
    t5 = t5.sort_values(by=["T5_signal"], ascending=False)
    t5.to_csv("t5_consensus_surfaceome.csv")
    df = t4.merge(t5, how="outer").fillna(0)
    df.to_csv("t4t5_consensus_surfaceome.csv", index=False)


@click.command()
@click.option("-i",
              "--input",
              help="Proteins table obtained after ratiometic classification")
def cli(input):
    "Command line interface"
    get_consensus_surfaceomes(input)


if __name__ == "__main__":
    cli()
