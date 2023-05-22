#!/usr/bin/env python3

"""
"""

import click
import json
import pandas as pd


def load_fubar(fubar):
    """
    """

    with open(fubar, "r") as fh:
        data = json.load(fh)
    headers = [h[0] for h in data["MLE"]["headers"]]
    content = data["MLE"]["content"]["0"]
    df = pd.DataFrame(content).drop([6, 7], axis=1)
    df.columns = headers

    return df


def compute_fubar_feaures(fubar):

    df = load_fubar(fubar)
    df["omega"] = df["beta"] / df["alpha"]
    # Mean dN, dS, Omega, fraction positive selected, fraction negative selected
    frac_positive = len(df)
    frac_negative = len(df)
    features = [df["beta"].mean(), df["alpha"].mean(),
                df["omega"].mean(), frac_positive, frac_negative]
    return features


def load_absrel(absrel):
    """
    """

    with open(absrel, "r") as fh:
        data = json.load(fh)

    return data


def parse_absrel_branches(absrel):
    """
    """

    data = load_absrel(absrel)
    selected_branches = []
    species = list(data["branch attributes"]["0"].keys())
    for specie in species:
        content = data["branch attributes"]["0"][specie]
        if content["Corrected P-value"] < 0.05:
            selected_branches.append(specie)

    return selected_branches


@click.command()
@click.option("-f",
              "--fubar",
              help="FUBAR results from Hyphy")
@click.option("-a",
              "--absrel",
              help="aBSREL results from Hyphy")
def cli(fubar, absrel):
    print(compute_fubar_feaures(fubar))
    print(parse_absrel_branches(absrel))


if __name__ == "__main__":
    cli()
