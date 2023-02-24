#!/usr/bin/env python3

import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from itertools import chain


def get_expressed_genes(df, cell):
    dff = df[cell]
    genes = dff[dff > 0.8].dropna().index

    return genes


def plot_non_expressed_proteins(data, gene2fbgn, proteome, neuron):

    # Load probability inference from RNAseq data
    df = pd.read_csv(data, sep="\t")
    # Change gene names to FlyBase identifiers
    ids = pd.read_csv(gene2fbgn, names=["gene", "fbgn"])
    ids = dict(zip(ids.gene, ids.fbgn))
    df.replace(ids, inplace=True)
    # Set identifiers as the new index
    df.set_index("Unnamed: 0", inplace=True)
    proteome = pd.read_csv(proteome)
    dff = df[df.index.isin(proteome.Protein)]
    non_expr_ref = dff[dff[neuron] < 0.2]
    sns.set(style="darkgrid", font_scale=1.6, rc={"lines.linewidth": 3})
    f, ax1 = plt.subplots(figsize=(8, 6))
    if neuron == "T4":
        color = "#f1a340"
    elif neuron == "T5":
        color = "#998ec3"
    sns.histplot(dff["T4"], ax=ax1, color=color)
    ax1.set_xlabel("Expression Probability")
    ax1.set_ylabel("# Proteins")
    plt.tight_layout()
    plt.savefig(f"{neuron}_dist.png")
    plt.close()

    return non_expr_ref, neuron


def get_connectome_counts(df, neuron):
    """
    """

    glia_types = ["Glia_Eg", "Glia_Mg", "Glia_Psg"]
    glia_genes = [get_expressed_genes(df, glia) for glia in glia_types]
    glia_genes = set(chain(*glia_genes))
    if neuron == "T4":
        input_types = ["Mi1", "Tm3", "Mi9", "Mi4"]
    elif neuron == "T5":
        input_types = ["Tm9", "Tm1", "Tm2"]
    input_genes = [get_expressed_genes(df, input) for input in input_types]
    input_genes = set(chain(*input_genes))
    output_types = ["LPC1", "LLPC1", "TmY5a", "LPLC1"]
    output_genes = [get_expressed_genes(df, output) for output in output_types]
    output_genes = set(chain(*output_genes))
    connectome = input_genes.union(output_genes)

    if neuron == "T4":
        colors = ("#fee391", "#f1a340")
    elif neuron == "T5":
        colors = ("#bcbddc", "#998ec3")
    venn2([glia_genes, connectome], ("Glia", "Connectome"), set_colors=colors, alpha=1)
    plt.tight_layout()
    plt.savefig(f"dummy{neuron}.png")
    plt.close()

# def plot_expression_data(df1, df2, output):
#     f, ax = plt.subplots()
#     sns.histplot(df1.iloc[:,0], stat="probability", binwidth=0.1, color="#f1a340", ax=ax)
#     sns.histplot(df2.iloc[:,0], stat="probability", binwidth=0.1, color="#998ec3", ax=ax)
#     ax.set_xlim(0, 3.5)
#     plt.savefig(f"{output}.png")


# Command line interface options
@click.command()
@click.option("-d", "--data",
              help="RNAseq data")
@click.option("-f", "--gene2fbgn",
              help="Mapping between gene name and FlyBase identifier")
@click.option("-t4", "--t4_proteome",
              help="T4 proteome list")
@click.option("-t5", "--t5_proteome",
              help="T5 proteome list")
def cli(data, gene2fbgn, t4_proteome, t5_proteome):
    """
    Command line interface
    """

    non_expr, neuron = plot_non_expressed_proteins(data, gene2fbgn, t4_proteome, "T4")
    get_connectome_counts(non_expr, neuron)
    non_expr, neuron = plot_non_expressed_proteins(data, gene2fbgn, t5_proteome, "T5")
    get_connectome_counts(non_expr, neuron)


if __name__ == '__main__':
    cli()
