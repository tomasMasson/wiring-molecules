#!/usr/bin/env python3

import click
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def read_t4t5_scrnaseq_data(folder, sample, t4_proteome, t5_proteome):
    adata = sc.read_10x_mtx(
            folder,
            var_names='gene_ids',
            prefix=f"{sample}_",
            cache=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    df = adata.to_df()
    # Neuronal expression (GFP and nSyb)
    df = df[(df["AGgn0000001"] > 1) & (df["FBgn0013342"] > 1)]
    # T4 neurons (TfAP-2 +)
    t4 = df[df["FBgn0261953"] > 0.8].mean()
    t4_proteome = pd.read_csv(t4_proteome)
    t4 = t4[t4.index.isin(t4_proteome.Protein)]
    expr_t4 = pd.DataFrame({"Expression": t4})
    expr_t4.loc[expr_t4["Expression"] <= 0.5, "Class"] = "Intermediate"
    expr_t4.loc[expr_t4["Expression"] <= 0.2, "Class"] = "Low"
    expr_t4.loc[expr_t4["Expression"] > 0.5, "Class"] = "High"
    expr_t4["Sample"] = f"{sample}"
    expr_t4["FBgn"] = expr_t4.index
    expr_t4.to_csv(f"{sample}_no_expr_t4.csv", index=False, header=False)
    # T5 neurons (TfAP-2 -)
    t5 = df[df["FBgn0261953"] <= 0.2].mean()
    t5_proteome = pd.read_csv(t5_proteome)
    t5 = t5[t5.index.isin(t5_proteome.Protein)]
    expr_t5 = pd.DataFrame({"Expression": t5})
    expr_t5.loc[expr_t5["Expression"] <= 0.5, "Class"] = "Intermediate"
    expr_t5.loc[expr_t5["Expression"] <= 0.2, "Class"] = "Low"
    expr_t5.loc[expr_t5["Expression"] > 0.5, "Class"] = "High"
    expr_t5["Sample"] = f"{sample}"
    expr_t5["FBgn"] = expr_t5.index
    expr_t5.to_csv(f"{sample}_no_expr_t5.csv", index=False, header=False)
    return expr_t4, expr_t5


def plot_expression_data(df1, df2, output):
    f, ax = plt.subplots()
    sns.histplot(df1.iloc[:,0], stat="probability", binwidth=0.1, color="#f1a340", ax=ax)
    sns.histplot(df2.iloc[:,0], stat="probability", binwidth=0.1, color="#998ec3", ax=ax)
    ax.set_xlim(0, 3.5)
    plt.savefig(f"{output}.png")


# Command line interface options
@click.command()
@click.option("-f", "--folder",
              help="Folder with the scRNAseq data")
@click.option("-s", "--sample",
              help="File name prefix")
@click.option("-t4", "--t4_proteome",
              help="T4 proteome list")
@click.option("-t5", "--t5_proteome",
              help="T5 proteome list")
def cli(folder, sample, t4_proteome, t5_proteome):
    """
    Command line interface
    """

    t4, t5 = read_t4t5_scrnaseq_data(folder, sample, t4_proteome, t5_proteome)
    plot_expression_data(t4, t5, sample)


if __name__ == '__main__':
    cli()

