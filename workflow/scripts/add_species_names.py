#!/usr/bin/env python3

import click
from Bio import SeqIO


def add_species_names(multifasta):
    """
    Adds species name to the protein identifier
    for multifastas coming from UniProt
    """

    # Load sequence data
    seqs = list(SeqIO.parse(multifasta, "fasta"))
    for s in seqs:
        specie = "_".join(s.description.split("OS=")[1].split()[:2])
        header = s.id + "|" + specie
        print(f">{header}\n{s.seq}")


@click.command()
@click.option("-i",
              "--input",
              help="Multifasta file to be processed")
def cli(input):
    "Command line interface"
    add_species_names(input)


if __name__ == "__main__":
    cli()
