#!/usr/bin/env python3

"Computes proteins members of each orthogroup at OrthoDB for a given taxon (e.g. Diptera)"

import click
import pandas as pd
from Bio import SeqIO


def get_protein_names(multifasta, output):
    """
    Creates a CSV file with the protein ID (Flybase) and its name
    """

    seqs = SeqIO.parse(multifasta, "fasta")
    with open(output, "w") as fh:
        for seq in seqs:
            description = seq.description
            protein_name = description.split(";")[3].split("name=")[1]
            fh.write(f"{seq.id},{protein_name}\n")



def filter_orthogroups(ogs, taxid):
    """
    Returns a filtered DataFrame from OrthoDB's orthogroups based on the TaxID
    """

    # Column names
    COLUMNS = ["OG", "Level", "OG_name"]
    # Load orthogroups table
    df = pd.read_csv(ogs,
                     sep="\t",
                     names=COLUMNS)
    # Filter orthogroups according to the Taxid
    fdf = df[df["Level"] == taxid]
    fdf.to_csv("diptera_orthogroups_intermediate.csv", index=False)

    return fdf


def map_ogs2genes(og2genes, ogs):
    """
    Returns a filtered gene list based on a set of orthogroups provided
    """

    # Genes table column names
    COLUMNS = ["OG", "Gene_ID"]
    # Initialize output dataframe
    df = pd.DataFrame()
    # Load genes table using smaller chunks
    input_iterator = pd.read_csv(og2genes,
                                 sep="\t",
                                 names=COLUMNS,
                                 chunksize=1000000)
    # Load orthogroups set
    ogs_set = ogs["OG"]
    # Iterate over genes table chunks
    for chunk in input_iterator:
        # Filter for genes inside orthogroups set
        chunk_flt = chunk[chunk["OG"].isin(ogs_set)]
        # Add genes to output dataframe
        df = pd.concat([df, chunk_flt], ignore_index=True)
    # Filter only orthogroups on the list
    fdf = df.merge(ogs, on=["OG"])

    return fdf


def add_gene_annotations(fgenes, annotations, species):
    """
    Add annotations and species information to the gene/orthogroup dataframe
    """

    # Read annotations table
    COLUMNS = ["Gene_ID",
               "OrthoDB_species"]
    # Initialize DataFrame for gene annotations
    df = pd.DataFrame()
    # Load genes table using smaller chunks
    ann_iterator = pd.read_csv(annotations,
                               sep="\t",
                               names=COLUMNS,
                               usecols=[0, 1],
                               chunksize=1000000)
    # Load list of filtered genes
    gene_list = list(fgenes["Gene_ID"])
    # Iterate over genes annotations chunks
    for chunk in ann_iterator:
        # Filter for annotations inside genes list
        chunk_flt = chunk[chunk["Gene_ID"].isin(gene_list)]
        # Add annotations to output dataframe
        df = pd.concat([df, chunk_flt], ignore_index=True)

    # Set column names
    COLUMNS = ["OrthoDB_species",
               "Scientific_name"]
    # Load species table
    species = pd.read_csv(species,
                          sep="\t",
                          names=COLUMNS,
                          usecols=[1, 2])

    # Merge genes and annotations DataFrames
    mdf = fgenes.merge(df, on=["Gene_ID"])
    # Merge annotated genes and species DataFrames
    mdf = mdf.merge(species, on=["OrthoDB_species"])
    # Write final table to CSV
    mdf.to_csv("diptera_orthogroups.csv", sep="\t")

    mdf_dmel = mdf[mdf["Scientific_name"] == "Drosophila melanogaster"]
    genes_map = pd.read_csv("diamond_search_dmel_proteome", sep="\t", usecols=[0, 1], names=["OG", "Gene_name"])
    mdf_dmel =mdf_dmel.set_index("Gene_ID").join(genes_map.set_index("Gene_ID"))
    names = get_protein_names("dmel_proteome", "dmel_proteins_names")
    mdf_dmel = pd.merge(mdf_dmel, names, on="Gene_name")

# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--orthogroups',
              help="Orthogroups table (e.g. 'odb10v1_OGs.tab')")
@click.option('--taxid',
              type=int,
              help="Taxon identifier used to filter OrthoDB orthologs table")
@click.option('--og2genes',
              help="Mapping from orthogroups to gene names")
@click.option('--genes',
              help="Table with genes cross database information")
@click.option('--species',
              help="Individual species information")
# CLI main function
def command_line_interface(orthogroups, taxid, og2genes, genes, species):
    """
    Computes proteins members of each orthogroup at OrthoDB for a given taxon
    """

    ogs = filter_orthogroups(orthogroups, taxid)
    ogs_filtered = map_ogs2genes(og2genes, ogs)
    add_gene_annotations(ogs_filtered, genes, species)


if __name__ == '__main__':
    command_line_interface()
