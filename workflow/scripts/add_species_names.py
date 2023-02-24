#!/usr/bin/env python3

import click
import pandas as pd
from pathlib import Path
from Bio import SeqIO


def add_species_names(multifasta, mapping):
    """
    Adds species name to the protein identifier
    for multifastas coming from UniProt
    """

    species = pd.read_csv(mapping, sep="\t")
    species_dic = dict(zip(species.Proteome, species.Id_Organism))
    p = Path(multifasta).stem.split("_")[0]
    specie = species_dic[p]
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
@click.option("-m",
              "--mapping",
              help="Name mapping between proteomes and species name")
def cli(input, mapping):
    "Command line interface"
    add_species_names(input, mapping)


if __name__ == "__main__":
    cli()
