#!/usr/bin/env python3

"""
Extract the members of each orthogroup stored at OrthoDB for a given taxon
"""

import click
import pandas as pd


def filter_orthogroups_table(ogs, taxid):
    """
    Returns a filtered DataFrame from OrthoDB orthogroups table (odb10v1_OGs)
    based on the TaxID provided
    """

    # Column names
    COLUMNS = ["OG", "Level", "OG_name"]
    # Load orthogroups table
    df = pd.read_csv(ogs,
                     sep="\t",
                     names=COLUMNS)
    # Filter orthogroups based on Taxid
    fdf = df[df["Level"] == taxid]

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
    # Load orthogroups list
    ogs_list = ogs["OG"]
    # Iterate over genes table chunks
    for chunk in input_iterator:
        # Keep genes present on orthogroups list
        chunk_flt = chunk[chunk["OG"].isin(ogs_list)]
        # Add genes to output dataframe
        df = pd.concat([df, chunk_flt], ignore_index=True)
    # Merge OG and OG2genes tables
    fdf = df.merge(ogs, on=["OG"])

    return fdf


def add_gene_annotations(fgenes, annotations, species):
    """
    Add gene annotations (odb10v1_genes) and species (odb10v1_species)
    information to the orthogroups dataframe
    """

    # Load genes table
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
    mdf.to_csv("orthogroups_table.csv", sep="\t", index=False)

    return mdf


def group_genes_by_orthogroup(df, mappings, output):
    """
    Compute the number and list of proteins belonging to each orthogroup
    """

    # Load mappings between OrthoDB and Flybase identifiers
    mappings = pd.read_csv(mappings,
                           sep="\t",
                           names=["orthodb", "flybase"])
    maps = dict(zip(mappings.orthodb, mappings.flybase))
    # Group genes by orthogroup
    ogs = df.groupby("OG")
    # dff = pd.DataFrame(columns=["Gene", "#_orthologs", "List_orthologs"])
    data = []
    COLUMNS = ["Gene_name", "#_orthologs", "List_orthologs"]
    for og in ogs:
        # Keep only orthogroups that contain D. melanogaster genes
        species = list(og[1]["Scientific_name"])
        if "Drosophila melanogaster" in species:
            # Store species name
            name = og[1][og[1].Scientific_name == "Drosophila melanogaster"].Gene_ID.values[0]
            # Store gene name
            genes = ",".join(list(og[1].Gene_ID))
            # Store number of genes in the orthogroup
            n_genes = len(og[1])
            data.append([name, n_genes, genes])

    dff = pd.DataFrame(data, columns=COLUMNS)
    dff = dff.replace(maps)
    dff.to_csv(output, sep="\t", index=False)
    print(dff)


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-o', '--orthogroups',
              help="Orthogroups table (e.g. 'odb10v1_OGs.tab')")
@click.option('-t', '--taxid',
              type=int,
              help="Taxon identifier used to filter OrthoDB orthologs table")
@click.option('-og2g', '--og2genes',
              help="Mapping from orthogroups to gene names")
@click.option('-g', '--genes',
              help="Table with genes cross database information")
@click.option('-s', '--species',
              help="Individual species information")
@click.option('-m', '--mappings',
              help="Mappigns from OrthoDB to Flybase identifiers")
@click.option('-o', '--output',
              help="File name for the output table")
# CLI main function
def cli(orthogroups, taxid, og2genes, genes, species, mappings, output):
    """
    Computes proteins members of each orthogroup at OrthoDB for a given taxon
    """

    ogs = filter_orthogroups_table(orthogroups, taxid)
    ogs_filt = map_ogs2genes(og2genes, ogs)
    df = add_gene_annotations(ogs_filt, genes, species)
    group_genes_by_orthogroup(df, mappings, output)


if __name__ == '__main__':
    cli()
