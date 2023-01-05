#!/usr/bin/env python3

import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def plot_venn_diagram(data1, data2):
    """
    Generate a two-way Venn diagram from T4 and T5 consensus protein lists
    """

    t4 = set(pd.read_csv(data1).Protein)
    t5 = set(pd.read_csv(data2).Protein)
    venn2([t4, t5], ("T4", "T5"), set_colors=("#f1a340", "#998ec3"), alpha=1)
    plt.tight_layout()
    plt.savefig("venn_diagram.png")

# Command line interface options
@click.command()
@click.option("-d1", "--data1",
              help="Consensus protein list for T4 samples")
@click.option("-d2", "--data2",
              help="Consensus protein list for T5 samples")
def cli(data1, data2):
    """
    Command line interface
    """

    plot_venn_diagram(data1, data2)


if __name__ == '__main__':
    cli()
