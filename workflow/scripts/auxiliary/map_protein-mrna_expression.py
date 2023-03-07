#!/usr/bin/env python3

"Module docstring"

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read mass-spec data
ms_df = pd.read_csv("pl_proteomics/T4T5_surface-ms_enriched_annotated.csv", sep="\t")

# Read RNA-seq data
mrna_df = pd.read_csv("pl_proteomics/GSE116969_dataTable7b.genes_x_cells_p_expression.modeled_genes.txt", sep="\t", index_col="Genes")

mrna_T4T5 = mrna_df[["T4", "T5", "T4.T5"]]
mrna_T4T5_exp = mrna_T4T5[mrna_T4T5 > 0.5].dropna()
genes_exp = list(mrna_T4T5_exp.index)

prot_exp = []
prot_not_exp = []

for prot in ms_df.SYMBOL:
    if prot in genes_exp:
        prot_exp.append(prot)
    else:
        prot_not_exp.append(prot)

pd.Series(prot_exp).to_csv('proteins_expressed.csv', index=False)
pd.Series(prot_not_exp).to_csv('proteins_not-expressed.csv', index=False)
