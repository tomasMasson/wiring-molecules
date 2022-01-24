#!/usr/bin/env python3

"Extract the evolutionary rate (dN/dS) from Hyphy FUBAR.json output file"

import click
import json
import numpy as np
import pandas as pd


def extract_fubar_evolrate(fubar, output):
    "Extract the evolutionary rate (dN/dS) from Hyphy FUBAR.json output file"

    # Open FUBAR.json file
    with open(fubar, 'r') as fh:
        data = fh.read()
    # Decode json
    data = json.loads(data)
    # Extract the computations from the data
    contents = data['MLE']['content']['0']
    # Set the header for the DataFrame
    headers = ['dS', 'dN']
    # Create a DataFrame
    df = pd.DataFrame(contents)
    # Filter out parameters that will not be used
    df = df.drop(columns=[2, 3, 4, 5, 6, 7])
    df.columns = headers
    df['Omega'] = np.divide(df['dN'], df['dS'])
    avg_ds = df['dS'].mean()
    avg_dn = df['dN'].mean()
    avg_omega = df['Omega'].mean()
    if output == 'average':
        print(avg_omega)
    elif output == 'profile':
        print(df.Omega.to_list())
    else:
        print('Please provide one of the following options as output: average or profile')



@click.command()
@click.option('-f', 'fubar')
@click.option('-o', 'output')
def cli(fubar, output):
    "Command line interface"
    extract_fubar_evolrate(fubar, output)


if __name__ == '__main__':
    cli()
