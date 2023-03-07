#!/usr/bin/env python3

import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def plot_venn_diagram(data):
    """
    Generate a two-way Venn diagram from T4 and T5 consensus protein lists
    """

    df = pd.read_csv(data)
    t4 = set(df[df["T4_signal"] != 0].Protein)
    t5 = set(df[df["T5_signal"] != 0].Protein)
    # t4 = set(pd.read_csv(data1, names=["Protein", "Signal"]).Protein)
    # t5 = set(pd.read_csv(data2, names=["Protein", "Signal"]).Protein)
    venn2([t4, t5], ("T4", "T5"), set_colors=("#f1a340", "#998ec3"), alpha=1)
    plt.tight_layout()
    plt.savefig("venn_diagram.svg")

# Command line interface options
@click.command()
@click.option("-d", "--data",
              help="Consensus protein list for T4 and T5 neurons")
def cli(data):
    """
    Command line interface
    """

    plot_venn_diagram(data)


if __name__ == '__main__':
    cli()
