#!/usr/bin/env python3

import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def plot_pca(data):
    """
    Generate a PCA plot for each sample based on the measured intensities
    """

    # Ingest data
    df = pd.read_csv(data, sep="\t").fillna(0)
    # Extract channel values and standarize them
    channels = StandardScaler().fit_transform(df.iloc[:, 2:].T)
    # Initialize PCA object
    pca = PCA(n_components=2)
    # Fit PCA model
    pca_components = pca.fit_transform(channels)
    # Load PCA into a DataFrame for plotting and add sample labels
    pca_df = pd.DataFrame(data=pca_components, columns=["PC1", "PC2"])
    pca_df["Sample"] = ["HRP_1", "H2O2_1", "HRP_2", "H2O2_2", "T4_1", "T4_2", "T5_1", "T5_2"]
    # Plot PCA scatterplot
    sns.set(style="darkgrid", font_scale=1.6, rc={"lines.linewidth": 3})
    f, ax1 = plt.subplots(figsize=(8, 6))
    sns.scatterplot(x="PC1", y="PC2", data=pca_df, hue="Sample", legend=False, s=180, palette=["#808080", "#808080", "#808080", "#808080", "#f1a340", "#f1a340", "#998ec3", "#998ec3"], ax=ax1)
    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")
    plt.tight_layout()
    plt.savefig("pca_analysis.png")

# Command line interface options
@click.command()
@click.option("-d", "--data",
              help="Raw data containing the spectral data")
def cli(data):
    """
    Command line interface
    """

    plot_pca(data)


if __name__ == '__main__':
    cli()
