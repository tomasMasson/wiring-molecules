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
    Returns a series with GO terms related with the specified localization
    """

    # Neuronal cell surface localization set
    surface = ["axons", "basal plasma membrane", "basolateral plasma membrane",
               "cell surface", "dendrite", "dendritic spine",
               "extracellular region", "extracellular space", "gap junction",
               "synapse", "plasma membrane", "neuron projection",
               "postsynaptic density", "postsynaptic membrane", "presynapse",
               "presynaptic membrane"]

    # Non surface localizations set
    non_surface = ["cytoplasm", "nucleus", "mitochondrial inner membrane",
                   "mitochondrial membrane", "mitochondrial matrix"]

    return (go_terms
            [go_terms["GO NAME"].isin(surface)]
            ["GENE PRODUCT ID"],
            go_terms
            [go_terms["GO NAME"].isin(non_surface)]
            ["GENE PRODUCT ID"]
            )


def load_go_terms(go_terms_table):
    """
    Ingests and filters the GO terms table to keep only ID and GENE PRODUCT names

    Returns a tuple with the lists for surface and non-surface proteins
    """

    df = pd.read_csv(go_terms_table,
                     sep="\t",
                     usecols=["GO NAME", "GENE PRODUCT ID"])
    return get_go_terms(df)


def map_identifiers(df, df2):

    ids = df.Protein_id.str.split(";")
    d = {}
    for recs in ids.iteritems():
        for rec in recs[1]:
            rec = rec.split('-')[0]
            mp = df2[df2.UniprotKB == rec].UniprotKB
            if len(mp) > 0:
                d[recs[0]] = list(mp)[0]
                break
    dff = df.rename(index=d)
    dff["Protein_id"] = dff.index
    print(dff.Protein_id)
    return dff


def dedup_protein_ids(df):
    """
    Transform a column with multiple proteins identifiers into a column with only one
    """

    return (df.assign(Protein_id=df.Protein_id
            .str.split(";", 1, expand=True)[0])
            )


def compute_tp_fp(data, map, sample, control, go_terms):
    """
    Calculates true positives and false positives for the sample and control
    """

    # Load GO terms data
    go_terms_tp, go_terms_fp = load_go_terms(go_terms)
    # Load TMT quantification and remove missing values
    return (data
            .dropna()
            .pipe(map_identifiers, map)
            # .pipe(dedup_protein_ids)
            .assign(signal_ratio=data[sample]-data[control],
                    tp=data.Protein_id.isin(go_terms_tp).astype("int"),
                    # tp=data.index.isin(go_terms_tp).astype("int"),
                    fp=data.Protein_id.isin(go_terms_fp).astype("int"),
                    # fp=data.index.isin(go_terms_fp).astype("int"),
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
                    tpr_fpr=(df.tp.cumsum() / df.tp.sum() - df.fp.cumsum() / df.fp.sum())
                    )
            )


def save_significant_proteins(df, sample, control):
    """
    """

    # print(df[df.signal_ratio > df.loc[df.tpr_fpr.idxmax(), "signal_ratio"]])
    (df
     [df.signal_ratio > df.loc[df.tpr_fpr.idxmax(), "signal_ratio"]]
     .to_csv(f"{sample}_{control}.csv", index=False)
     )

    return df


def plot_ratiometric_analysis(df, sample, control):
    """
    pass
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
    # renderer = hv.plotting.mpl.MPLRenderer.instance(dpi=400)
    # renderer.save(p, f'{sample}_{control}', fmt='svg')
    # p_plot = renderer.get_plot(p)
    hv.save(p, f'{sample}_{control}.png')


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


# Command line interface options
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-d", "--data",
              help="TMT quantification data")
@click.option("-m", "--map",
              help="Identifiers mapping")
@click.option("-s", "--sample",
              help='Sample to run the analysis ["t4_1", "t4_2", "t5_1", "t5_2"]')
@click.option("-c", "--control",
              help='Control to normalize the analysis ["hrp_1", "h2o2_1", "hrp_2", "h2o2_2"]')
@click.option("-g", "--go_terms",
              help="GO terms table")
# CLI main function
def cli(data, map, sample, control, go_terms):
    """
    Command line interface to set ratiometric analysis and plot data
    """

    data = pd.read_csv(data, sep="\t")
    map = pd.read_csv(map)
    return (data
            .pipe(compute_tp_fp, map, sample, control, go_terms)
            .pipe(compute_tpr_fpr)
            .pipe(save_significant_proteins, sample, control)
            .pipe(plot_ratiometric_analysis, sample, control)
            )


if __name__ == '__main__':
    cli()