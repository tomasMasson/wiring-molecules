#!/usr/bin/env python3

import click
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def plot_go_terms_enrichment(data, go_class):
    """
    """

    # Ingest data
    df = pd.read_csv(data)
    df = df[df["GO-Class"] == go_class]
    df["-log10(FDR)"] = -np.log10(df["FDR"])
    t4 = df[df["Neuron"] == "T4"]
    t5 = df[df["Neuron"] == "T5"]
    t5["-log10(FDR)"] = t5["-log10(FDR)"].mul(-1)
    dff = pd.concat([t4, t5])
    # Plot GO terms barplot
    sns.set(style="darkgrid", font_scale=3,
            rc={"lines.linewidth": 2})
    f, ax1 = plt.subplots(figsize=(14, 10))
    p = sns.barplot(x="-log10(FDR)", y="Description", data=dff, palette=["#f1a340", "#998ec3"], hue="Neuron", dodge=False, ax=ax1)
    p.legend_.remove()
    p.set_ylabel("")
    if go_class == "BP":
        ax1.set_title("Biological Process",
                      weight="bold")
    if go_class == "CC":
        ax1.set_title("Cellular Component",
                      weight="bold")
    plt.tight_layout()
    plt.savefig(f"stringdb_{go_class}_analysis.svg")


# Command line interface options
@click.command()
@click.option("-d", "--data",
              help="Raw data containing StringDB GO terms annotations")
@click.option("-c", "--go_class",
              help="GO class to be plotted")
def cli(data, go_class):
    """
    Command line interface
    """

    plot_go_terms_enrichment(data, go_class)


if __name__ == '__main__':
    cli()
