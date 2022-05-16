#!/usr/bin/env python3

"""
Process TMT mass spectrometry abundance values
and plot the ROC curves for the samples.
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import holoviews as hv
hv.extension("bokeh")
from scipy.stats import ranksums


def compute_signal_ratio(sample, control):
    """
    Normalize sample intensity with respect to some
    arbitrary control. Signal ratio is expressed as
    log10.
    """

    s = pd.to_numeric(sample)
    c = pd.to_numeric(control)
    ratio = s.subtract(c)

    return ratio


def rank_dataframe(df, key):
    """
    Sorts a DataFrame in descending order based on key
    parameter
    """

    dfr = df.sort_values(by=[key], ascending=False)
    return dfr


def plot_ratiometric_analysis(df, s_name, c_name):
    """
    pass
    """

    # Set color palette
    TP_COLOR = "#fdb462"
    FP_COLOR = "#000000"
    # Define reference line for ROC curve
    line = pd.DataFrame({"x": np.linspace(0, 1, 10),
                         "y": np.linspace(0, 1, 10)})
    # Plot ROC curve
    roc = hv.Curve((df.FPR, df.TPR),
                   ("FPR"),
                   ("TPR"),
                   group="ROC Curve",
                   label=f"{s_name}")
    roc.opts(color=TP_COLOR)
    roc_random = hv.Curve((line.x, line.y),
                          ("FPR"),
                          ("TPR"),
                          group="ROC Curve",
                          label="Random classifier")
    roc_random.opts(color="grey",
                    line_dash="dotted")

    # Get true positive data points
    tp = df[df["TP"] == 1]
    frequencies, edges = np.histogram(tp.Signal_ratio, bins=20, range=(-1.5, 1.5), density=True)
    tp_dist = hv.Histogram((edges, frequencies))
    tp_dist.opts(color=TP_COLOR,
                 alpha=0.5,
                 xlabel=f"log10({s_name}/{c_name})",
                 ylabel="Density",
                 xlim=(-1.5, 1.5))
    # Get false positive data points
    fp = df[df["FP"] == 1]
    frequencies, edges = np.histogram(fp.Signal_ratio, bins=20, range=(-1.5, 1.5), density=True)
    fp_dist = hv.Histogram((edges, frequencies))
    fp_dist.opts(color=FP_COLOR,
                 alpha=0.5,
                 xlabel=f"log10({s_name}/{c_name})",
                 ylabel="Density",
                 xlim=(-1.5, 1.5))

    tpr_fpr = hv.Curve((df.Signal_ratio, df["TPR-FPR"]))
    tpr_fpr.opts(color=TP_COLOR,
                 xlabel=f"log10({s_name}/{c_name})",
                 ylabel="TPR-FPR")

    r = roc * roc_random
    r.opts(legend_position="bottom_right")
    p = r + tp_dist * fp_dist + tpr_fpr
    hv.save(p, "plot.html", backend="bokeh")

    # fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    # # Compute Area Under the Curve (AUC) value
    # auc = np.trapz(df.TPR, df.FPR)
    # # Calculate Wilcoxon rank-sum significance value
    # pvalue = np.round(ranksums(df.TPR, df.FPR, alternative="greater")[1], 3)
    # # Plot receiver operating characteristic curve
    # line = pd.DataFrame({"x": np.linspace(0, 1, 10),
    #                      "y": np.linspace(0, 1, 10)})
    # sns.lineplot(data=df, x="FPR", y="TPR", ax=ax1)
    # sns.lineplot(data=line, x="x", y="y", ax=ax1)
    # ax1.text(0.1, 1.0, pvalue)
    # # Plot True and False positives distributions
    # dff = df[df["TP"] == 1]
    # sns.histplot(data=dff.Signal_ratio,
    #              stat="probability",
    #              color="#1f78b4",
    #              bins=30,
    #              ax=ax2)
    # dff = df[df["FP"] == 1]
    # sns.histplot(data=dff.Signal_ratio,
    #              stat="probability",
    #              color="#b2df8a",
    #              bins=30,
    #              ax=ax2)
    # # Plot TRP-FPR vs intensity ratio curve
    # sns.lineplot(data=df, x="Signal_ratio", y="TPR-FPR", ax=ax3)
    # # plt.show()
    # plt.savefig(f"{sample_name}.png")
    # plt.close()


def compute_ratiometric_analysis(data, sample, control, go_terms):
    """
    Calculates the True Positive and False Positive rates
    (TPR and FPR, respectively) for the sample and control.
    """

    # Load TMT quantification
    df = pd.read_csv(data, sep="\t")
    # Discard missing values
    df = df.dropna()
    # Assign an unique IDs for each protein
    df["Protein_unique_id"] = df["Protein_id"].str.split(";", 1, expand=True)[0]

    # Load GO terms annotations
    go_terms = pd.read_csv(go_terms, sep="\t")
    # True positive terms set
    tp_list = ["axons", "basal plasma membrane", "basolateral plasma membrane", "cell surface", "dendrite", "dendritic spine", "extracellular region", "extracellular space", "gap junction", "synapse", "plasma membrane", "neuron projection", "postsynaptic density", "postsynaptic membrane", "presynapse", "presynaptic membrane"]
    go_terms_tp = go_terms[go_terms["GO NAME"].isin(tp_list)]
    # False positive terms set
    # fp_list = ["cytoplasm", "nucleus", "mitochondrial inner membrane", "mitochondrial membrane", "mitochondrial matrix"]
    fp_list = ["cytoplasm", "nucleus", "mitochondria", "mitochondrial matrix"]
    go_terms_fp = go_terms[go_terms["GO NAME"].isin(fp_list)]
    # Compute TP and FP in sample data
    df["TP"] = df["Protein_unique_id"].isin(go_terms_tp["GENE PRODUCT ID"]).astype("int")
    df["FP"] = df["Protein_unique_id"].isin(go_terms_fp["GENE PRODUCT ID"]).astype("int")

    # Compute signal ratio
    df["Signal_ratio"] = compute_signal_ratio(df[sample], df[control])
    # Rank DataFrame based on intensity ratios
    df = rank_dataframe(df, "Signal_ratio")
    # Compute TPR and FPR counts and differences
    df["TPR"] = df["TP"].cumsum() / df["TP"].sum()
    df["FPR"] = df["FP"].cumsum() / df["FP"].sum()
    df["TPR-FPR"] = df["TPR"] - df["FPR"]
    # Keep only relevant columns
    dff = df.loc[:, ["Protein_unique_id", "Gene_name", "TP", "FP", "TPR", "FPR", "TPR-FPR", "Signal_ratio"]]
    # Get signal ratio cutoff from the maximal difference between TPR and FPR
    idxmax = df["TPR-FPR"].idxmax()
    tpr_fpr_max = dff.loc[idxmax, "TPR-FPR"]
    sample_name = f"{sample}_{control}"
    plot_ratiometric_analysis(dff, sample, control)
    dff["Significant"] = df["Signal_ratio"] > tpr_fpr_max
    dff.loc[:, "Sample"] = sample_name
    dff.to_csv(f"{sample_name}.csv", index=False)

    return dff


# for sample in ["t4_1", "t4_2", "t5_1", "t5_2"]:
for sample in ["t4_1"]:
    # for control in ["hrp_1", "h2o2_1", "hrp_2", "h2o2_2"]:
    for control in ["hrp_1"]:
        print(sample, control)
        DATA = "resources/t4t5_mass-spec_dataset.csv"
        SAMPLE = sample
        CONTROL = control
        GO_TERMS = "resources/QuickGO-annotations-1652225849796-20220510.tsv"
        compute_ratiometric_analysis(DATA, SAMPLE, CONTROL, GO_TERMS)
