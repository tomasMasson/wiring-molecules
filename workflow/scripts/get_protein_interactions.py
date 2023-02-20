#!/usr/bin/env python3

import click
import pandas as pd


def get_protein_interactions(proteins, mappings, gene_names, interactome, output):
    """
    Starting from a reference interactome (in StringDB format),
    the script will filter only the interactions involving two proteins
    present on the input list.
    Protein and Gene unique ids will be mapped into their gene name according
    to FlyBase.

    Invocation example:

    ./get_protein_interactions.py --proteins {protein_list.csv} --mappings {fbgn_fbtr_fbpp_fb_2023_01.tsv} --gene_names {fbgn_annotation_ID.tsv} --interactome {7227.protein.links.v11.5.txt} --output {proteins_interactions.tsv}

    """

    # Load protein list (as FBgn identifiers)
    prots = pd.read_csv(proteins, names=["Protein"])

    # Load mappings from FBgn -> FBpp
    # Since the protein list is in gene format (FBgn),
    # we need to convert them to protein entries (FBpp)
    gene2prot = pd.read_csv(mappings,
                            skiprows=5,
                            names=["Gene",
                                   "Transcript",
                                   "Protein"],
                            sep="\t")
    # Keep only the proteins that are present on the initial dataset
    gene2prot = gene2prot[gene2prot["Gene"].isin(prots["Protein"])]

    # Load reference interactome
    df = pd.read_csv(interactome, sep=" ")
    # Rename proteins so they match FlyBase style
    df["protein1"] = df["protein1"].str.split(".", expand=True)[1]
    df["protein2"] = df["protein2"].str.split(".", expand=True)[1]

    # Filter interaction that involve only proteins present on the initial list
    ints = df[(df["protein1"].isin(gene2prot["Protein"])) & (df["protein2"].isin(gene2prot["Protein"]))]
    # Filter high-confidence interactions (StringDB score >= 700)
    ints = ints[ints.combined_score >= 700]

    # Add interaction type column
    ints["Type"] = "pp"
    ints = ints[["protein1", "Type", "protein2"]]
    # Add aesthetic interactions for rendering the protein network on Cytoscape
    prots.columns = ["protein1"]
    prots["Type"] = "pn"
    prots["protein2"] = "T4/T5"
    ints = pd.concat([ints, prots])

    # Replace FBpp identifiers with the gene names
    ints.replace(dict(zip(gene2prot.Protein, gene2prot.Gene)),
                 inplace=True)
    # Load gene names file
    genes = pd.read_csv(gene_names, sep="\t",
                        skiprows=5,
                        names=["Gene", "Org", "FBgn"],
                        usecols=[0, 1, 2])
    # Initialize a dictionary to replace FBgn -> Gene
    gene_dic = dict(zip(genes["FBgn"], genes["Gene"]))
    ints.replace(gene_dic, inplace=True)

    # Save interactions into output
    ints.to_csv(output, sep="\t",
                index=False, header=True)


# Command line interface options
@click.command()
@click.option("-p", "--proteins",
              help="Proteins list to filter interactions")
@click.option("-m", "--mappings",
              help="Gene to protein (FBgn -> FBpp) mappings as downloaded from FlyBase")
@click.option("-g", "--gene_names",
              help="Gene names to replace FBgn identifiers")
@click.option("-i", "--interactome",
              help="Reference interactome (from StringDB) to get protein interactions")
@click.option("-o", "--output",
              help="Output file name")
def cli(proteins, mappings, gene_names, interactome, output):
    """
    Command line interface
    """

    get_protein_interactions(proteins, mappings, gene_names, interactome, output)


if __name__ == "__main__":
    cli()
