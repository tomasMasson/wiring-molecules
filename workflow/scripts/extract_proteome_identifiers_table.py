#!/usr/bin/env python3

import click
import pandas as pd
from Bio import SeqIO


def parse_protein_identifiers(seq):
    """
    Takes a BioPython seqrecord and parse the Flybase and UniprotKB
    identifiers from the fasta description
    """

    # Load fasta header description
    desc = seq.description
    # Parse Uniprot id
    uniprot = desc.split(';')[5].split('dbxref=')[-1].split('UniProt')[-1].split(':')[1].split(',')[0]
    # Parse protein name
    gene_name = desc.split(';')[3].split('name=')[-1]
    # Parse FlyBase protein identifier
    fbpp = desc.split(';')[2].split('ID=')[-1]
    # Parse FlyBase gene identifier
    fbgn = desc.split(';')[4].split('parent=')[-1].split(',')[0]
    # Parse FlyBase transcript identifier
    fbtr = desc.split(';')[4].split('parent=')[-1].split(',')[1]

    return uniprot, gene_name, fbgn, fbtr, fbpp


def extract_protein_identifiers(multifasta, output):
    """
    Output a table with FlyBase and UniprotKB identifiers
    from a proteome
    """

    # Load proteome data
    seqs = SeqIO.parse(multifasta, "fasta")
    # Initialize lists to store identifiers
    fbgn = []
    fbtr = []
    fbpp = []
    gene_name = []
    uniprot = []

    # Extract data from each sequence
    for seq in seqs:
        unip, name, gn, tr, pp = parse_protein_identifiers(seq)
        # Pool data from each sequence
        uniprot.append(unip)
        gene_name.append(name)
        fbgn.append(gn)
        fbtr.append(tr)
        fbpp.append(pp)

    # Merge all data into a DataFrame
    data = {"UniprotKB": uniprot,
            "Gene": gene_name,
            "FBgn": fbgn,
            "FBtr": fbtr,
            "FBpp": fbpp
            }
    df = pd.DataFrame(data=data)
    # Save into tabular format
    df.to_csv(output, index=False)


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
# Command line interface options
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-s", "--seq",
              help="Sequence file used to extract proteins identifiers")
@click.option("-o", "--output",
              help="Output file")
# CLI main function
def cli(seq, output):
    """
    Command line interface
    """

    extract_protein_identifiers(seq, output)


if __name__ == "__main__":
    cli()
