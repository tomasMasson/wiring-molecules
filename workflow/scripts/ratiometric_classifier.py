#!/usr/bin/env python3

"""
Process TMT mass spectrometry abundance values
and plot ROC curves for the samples.
"""

import click
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ranksums


def get_subcellular_localization_references(go_terms, flyxcdb):
    """
    Returns the true and negative localization proteins sets
    """

    df = pd.read_csv(go_terms, sep="\t",
                     skiprows=5,
                     usecols=[0, 1, 2, 3, 4],
                     names=["DB", "DB_Object_ID",
                            "DB_Object_symbol",
                            "Qualifier", "GO_ID"])

    pm = set(df[df.GO_ID == "GO:0005886"].DB_Object_ID)
    ext = set(df[df.GO_ID == "GO:0005576"].DB_Object_ID)
    nuc = set(df[df.GO_ID == "GO:0005634"].DB_Object_ID)
    nucl = set(df[df.GO_ID == "GO:0005730"].DB_Object_ID)
    cyt = set(df[df.GO_ID == "GO:0005737"].DB_Object_ID)
    mit = set(df[df.GO_ID == "GO:0005739"].DB_Object_ID)
    flyxcdb = set(pd.read_csv(flyxcdb).GeneID)

    tp = pm.union(ext).union(flyxcdb)
    fp = nuc.union(nucl).union(cyt).union(mit)
    duplicates = tp.intersection(fp)
    fp = fp.difference(duplicates)
    tp = pd.DataFrame({"protein": list(tp),
                       "annotation": ["TP"] * len(tp)})
    fp = pd.DataFrame({"protein": list(fp),
                       "annotation": ["FP"] * len(fp)})
    go_terms = tp.append(fp)

    return (go_terms
            [go_terms["annotation"] == "TP"]
            ["protein"],
            go_terms
            [go_terms["annotation"] == "FP"]
            ["protein"]
            )


def map_identifiers(data, mappings):

    # Initialize a dict to store mappings
    map_dic = {}
    # Scan all the protein ids on the table
    for recs in data.Protein_id.iteritems():
        for rec in recs[1].split(";"):
            rec = rec.split("-")[0]
            # Map to FlyBase
            mp = mappings[mappings.uniprot == rec].flybase
            # Save positive results to the dict
            if len(mp) > 0:
                map_dic[recs[0]] = list(mp)[0]
                break

    # Output results as a DataFrame with a new "Protein_id" column
    return (data
            .drop("Protein_id", axis=1)
            .rename(index=map_dic)
            .reset_index()
            .rename(columns={"index": "Protein_id"})
            )


def compute_tp_fp(data, mappings, go_terms, flybase, sample, control):
    """
    Calculates true positives and false positives for the sample and control
    """

    # Load TMT quantification and map UniProt accession to FlyBase
    data = map_identifiers(data.dropna(), mappings)
    # Load GO terms data
    tp, fp = get_subcellular_localization_references(go_terms, flybase)
    # Compute TP, FP and signal ratios
    return (data
            .dropna()
            .assign(signal_ratio=data[sample]-data[control],
                    fp=data.Protein_id.isin(fp).astype("int"),
                    tp=data.Protein_id.isin(tp).astype("int"),
                    )
            .sort_values(by=["signal_ratio"], ascending=False)
            )


def compute_tpr_fpr(df):
    """
    Calculates true positives and false positives for the sample and control
    """

    return (df
            .assign(tpr=df.tp.cumsum() / df.tp.sum(),
                    fpr=df.fp.cumsum() / df.fp.sum(),
                    tpr_fpr=(df.tp.cumsum() / df.tp.sum() - df.fp.cumsum() / df.fp.sum()),
                    fdr=(df.fp.cumsum() / (df.tp.sum() + df.fp.sum()))
                    )
            )


def save_significant_proteins(df, sample, control):
    """
    Write the proteins that are above the maximum TPR-FPR threshold
    """

    (df
     [df.fdr < 0.1]
     .drop_duplicates(subset=["Protein_id"])
     .loc[:, ["Protein_id"]]
     .assign(Sample=f"{sample}")
     .assign(Control=f"{control}")
     .to_csv(f"{sample}_{control}.csv", index=False)
     )

    return df


def plot_ratiometric_analysis(df, sample, control):

    """
    Save diagnostic plot for the ratiometric analysis
    """

    # Set color palette
    TP_COLOR = "#f1a340"
    FP_COLOR = "#998ec3"
    # Compute pvalue and AUC score
    pvalue = str(np.round(ranksums(df.tpr, df.fpr, alternative="greater")[1], 2))
    auc = str(np.round(np.trapz(df.tpr, df.fpr), 4))
    # Define curve for random classifier
    line = pd.DataFrame({"x": np.linspace(df.fpr.min(), df.fpr.max(), 10),
                         "y": np.linspace(df.tpr.min(), df.tpr.max(), 10)})

    # Increase font size
    sns.set(style="darkgrid", font_scale=1.6, rc={"lines.linewidth": 3})
    # Define panels
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16,5))
    # Plot ROC curve for sample
    sns.lineplot(x=df.fpr, y=df.tpr,
                 color=TP_COLOR,
                 ax=ax1)
    sns.lineplot(x=line.x, y=line.y,
                 linestyle="--",
                 color=FP_COLOR,
                 ax=ax1)
    ax1.set_xlabel("False Positive Rate (FPR)")
    ax1.set_ylabel("True Postive Rate (TPR)")
    ax1.text(0.04, 0.94, f"AUC={auc}")
    # Plot TP and FP distributions
    tp =df[df.tp == 1].signal_ratio
    fp =df[df.fp == 1].signal_ratio
    sns.histplot(tp, bins=20,
                 color="#f1a340",
                 binwidth=0.1,
                 stat="probability",
                 ax=ax2)
    sns.histplot(fp, bins=20,
                 color="#998ec3",
                 binwidth=0.1,
                 stat="probability",
                 ax=ax2)
    ax2.set_xlabel("Signal Intensity Ratio")
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(0, 0.34)
    # Plot TPR-FPR distribution
    sns.lineplot(x=df.signal_ratio, y=df.tpr_fpr,
                 color="#f1a340",
                 ax=ax3)
    ax3.set_xlabel("Signal Intensity Ratio")
    ax3.set_ylabel("TPR - FPR")
    ax3.set_xlim(-2, 2)
    ax3.set_ylim(-0.1, 0.6)
    # Fix legend overlap
    plt.tight_layout(h_pad=2)
    fig.savefig(f"{sample}_{control}.png")


def run_analysis(data, mappings, annotations, flyxcdb, label):
    """
    Ratiometric analysis and plot data
    """

    # Load the data, mappings and localization datasets
    data = pd.read_csv(data, sep="\t")
    mappings = pd.read_csv(mappings,
                           sep="\t",
                           names=["uniprot", "flybase"])
    sample = "_".join(label.split("_")[0:2])
    control = "_".join(label.split("_")[2:])

    # Perform ratiometric analysis
    return (data
            .pipe(compute_tp_fp, mappings, annotations, flyxcdb, sample, control)
            .pipe(compute_tpr_fpr)
            .pipe(save_significant_proteins, sample, control)
            .pipe(plot_ratiometric_analysis, sample, control)
            )


# Command line interface options
@click.command()
@click.option("-d", "--data",
              help="Raw data containing the spectral data")
@click.option("-m", "--mappings",
              help="List detailing the mapping between UniProt IDs (present on the original table) and FlyBase IDs, employed throughout the analysis")
@click.option("-g", "--go_terms",
              help="Drosophila melanogaster proteomeGO terms mapping used to define True Positives (TP) and False Positives (FP) protein list based on subcellular localization")
@click.option("-f", "--flyxcdb",
              help="Drosophila melanogaster FlyXCDB is used to expand the set of True Positives (TP) protein list based on subcellular localization")
@click.option("-l", "--label",
              help="Label containing the sample and control names to be used in the analysis")
def cli(data, mappings, go_terms, flyxcdb, label):
    """
    Command line interface
    """

    run_analysis(data, mappings, go_terms, flyxcdb, label)


if __name__ == '__main__':
    cli()
