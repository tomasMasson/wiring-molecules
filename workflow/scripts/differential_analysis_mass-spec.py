#!/usr/bin/env python3

"""
Process TMT mass spectrometry abundance values
and plot ROC curves for the samples.
"""

import click
import numpy as np
import pandas as pd
import holoviews as hv
from scipy.stats import ranksums
hv.extension("bokeh")
hv.renderer('matplotlib')


def get_go_terms(go_terms):
    """
    Returns the true and negative localization proteins sets
    """

    return (go_terms
            [go_terms["annotation"] == "TP"]
            ["protein"],
            go_terms
            [go_terms["annotation"] == "FP"]
            ["protein"]
            )


def map_identifiers(df, df2):

    # Initialize a dict to store mappings
    map_dic = {}
    # Scan all the protein ids on the table
    for recs in df.Protein_id.iteritems():
        for rec in recs[1].split(";"):
            # Map to FlyBase
            mp = df2[df2.uniprot == rec].flybase
            # Save positive results to the dict
            if len(mp) > 0:
                map_dic[recs[0]] = list(mp)[0]
                break

    # Output results as a DataFrame with a new "Protein_id" column
    return (df
            .drop("Protein_id", axis=1)
            .rename(index=map_dic)
            .reset_index()
            .rename(columns={"index": "Protein_id"})
            )


def compute_tp_fp(data, mappings, localization, sample, control):
    """
    Calculates true positives and false positives for the sample and control
    """

    # Load TMT quantification and map UniProt accession to FlyBase
    data = map_identifiers(data.dropna(), mappings)
    # Load GO terms data
    go_terms_tp, go_terms_fp = get_go_terms(localization)
    # Compute TP, FP and signal ratios
    return (data
            .dropna()
            .assign(signal_ratio=data[sample]-data[control],
                    fp=data.Protein_id.isin(go_terms_fp).astype("int"),
                    tp=data.Protein_id.isin(go_terms_tp).astype("int"),
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
     .to_csv(f"{sample}_{control}.csv", index=False)
     )

    return df


def plot_ratiometric_analysis(df, sample, control):
    """
    Save diagnostic plot for the ratiometric analysis
    """

    # Set color palette
    TP_COLOR = "#fdb462"
    FP_COLOR = "#000000"
    # Compute pvalue and AUC score
    pvalue = str(np.round(ranksums(df.tpr, df.fpr, alternative="greater")[1], 2))
    auc = str(np.round(np.trapz(df.tpr, df.fpr), 4))
    # Define curve for random classifier
    line = pd.DataFrame({"x": np.linspace(0, 1, 10),
                         "y": np.linspace(0, 1, 10)})
    # Plot ROC curve for sample
    roc = hv.Curve((df.fpr, df.tpr),
                   ("FPR"),
                   ("TPR"),
                   group="ROC Curve",
                   label=f"{sample}")
    roc.opts(color=TP_COLOR)
    # Plot ROC curve for random classifier
    roc_random = hv.Curve((line.x, line.y),
                          ("FPR"),
                          ("TPR"),
                          group="ROC Curve",
                          label="Random classifier")
    roc_random.opts(color="grey",
                    line_dash="dotted")

    # Get true positive data points
    tp = df[df.tp == 1]
    frequencies, edges = np.histogram(tp.signal_ratio, bins=20, range=(-1.5, 1.5), density=True)
    tp_dist = hv.Histogram((edges, frequencies))
    tp_dist.opts(color=TP_COLOR,
                 alpha=0.5,
                 xlabel=f"log10({sample}/{control})",
                 ylabel="Density",
                 xlim=(-1.5, 1.5))
    # Get false positive data points
    fp = df[df.fp == 1]
    frequencies, edges = np.histogram(fp.signal_ratio, bins=20, range=(-1.5, 1.5), density=True)
    fp_dist = hv.Histogram((edges, frequencies))
    fp_dist.opts(color=FP_COLOR,
                 alpha=0.5,
                 xlabel=f"log10({sample}/{control})",
                 ylabel="Density",
                 xlim=(-1.5, 1.5))

    tpr_fpr = hv.Curve((df.signal_ratio, df["tpr_fpr"]))
    tpr_fpr.opts(color=TP_COLOR,
                 xlabel=f"log10({sample}/{control})",
                 ylabel="TPR-FPR")

    r = roc * roc_random * hv.Text(0.3, 1.0, f"AUC={auc}")
    r.opts(legend_position="bottom_right")
    p = r + tp_dist * fp_dist + tpr_fpr
    hv.save(p, f'{sample}_{control}.png')


def run_analysis(data, mappings, annotations, label):
    """
    Ratiometric analysis and plot data
    """

    # Load the data, mappings and localization datasets
    data = pd.read_csv(data, sep="\t")
    mappings = pd.read_csv(mappings,
                           sep="\t",
                           names=["uniprot", "flybase"])
    localizations = pd.read_csv(annotations)
    sample = "_".join(label.split("_")[0:2])
    control = "_".join(label.split("_")[2:])

    # Perform ratiometric analysis
    return (data
            .pipe(compute_tp_fp, mappings, localizations, sample, control)
            .pipe(compute_tpr_fpr)
            .pipe(save_significant_proteins, sample, control)
            .pipe(plot_ratiometric_analysis, sample, control)
            )


# Command line interface options
@click.command()
@click.option("-d", "--data",
              help="TMT quantification data")
@click.option("-m", "--mappings",
              help="Identifiers mapping")
@click.option("-a", "--annotations",
              help="Reference list with True Positives (TP) and False Positives (FP) for protein localizations")
@click.option("-l", "--label",
              help="Label containing the sample and control names to be used in the analysis")
def cli(data, mappings, annotations, label):
    """
    Command line interface
    """

    run_analysis(data, mappings, annotations, label)


if __name__ == '__main__':
    cli()
