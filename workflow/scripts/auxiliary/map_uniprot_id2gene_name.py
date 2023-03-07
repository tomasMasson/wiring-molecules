#!/usr/bin/env python3

"Module docstring"

import pandas as pd

uniprot_ids = pd.read_csv("data/fbgn_NAseq_Uniprot_fb_2021_04.tsv",
                         sep="\t",
                         names=["gene_symbol", "organism_abbreviation",
                                "primary_FBgn", "nucleotide_accession",
                                "na_based_protein_accession", "UniprotKB",
                                "EntrezGene_ID", "RefSeq_transcripts",
                                "RefSeq_proteins"],
                         skiprows=5,
                         skipfooter=1)
uniprot_ids = uniprot_ids.loc[:, ["gene_symbol", "UniprotKB"]]


gene_list = pd.read_csv("pl_proteomics/t4t5_enriched_proteins.csv",
                        names=["Protein"])

# Map uniprot identifiers to FlyBase gene names
genes_names = uniprot_ids.dropna().loc[uniprot_ids.dropna().UniprotKB.isin(gene_list.Protein), :]
genes_names.to_csv("t4t5_enriched_genes.csv", index=False)
