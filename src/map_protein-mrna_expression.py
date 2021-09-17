#!/usr/bin/env python3

"Module docstring"

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ms_df = pd.read_csv("pl_proteomics/t4t5_enriched_genes.csv",
                    names=["gene_symbol", "UniprotKB"])
# Pick upregulated genes in the MS dataset
upregulated_prot = ms_df["gene_symbol"]

mrna_df = pd.read_csv("pl_proteomics/GSE116969_dataTable7b.genes_x_cells_p_expression.modeled_genes.txt", sep="\t", index_col="Genes")
# Filter only T4 and T5 neurons
mrna_df = mrna_df.loc[:, ["T4", "T5"]]

mrna_ms_df = mrna_df.loc[mrna_df.index.isin(upregulated_prot), :]

t4 = list(mrna_ms_df[mrna_ms_df.T4 > 0.5].index)

sns.displot(data=mrna_ms_df.T4)
sns.displot(data=mrna_ms_df.T5)
plt.show()
