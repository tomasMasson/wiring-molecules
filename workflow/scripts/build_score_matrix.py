#!/usr/bin/env python3

"Module docstring"

from pathlib import Path
from Bio import SeqIO 
import pandas as pd

# Collect the rows names
seqs = SeqIO.parse("phylogenetic_profiles/UP000000803_7227.fasta", "fasta")
rows = [seq.id for seq in seqs]

# Collect the columns names
p = Path("phylogenetic_profiles/proteomes/")
columns = [child.name.split(".")[0] for child in p.iterdir()]

# Initialize an empty DataFrame
df = pd.DataFrame()

p = Path("phylogenetic_profiles/scores/")
for child in p.iterdir():
    # Store proteome name
    proteome = child.name.split(".")[0]
    names = ["qseqid", "sseqid", "pident", "length",
             "mismatch", "gapopen", "qstart", "qend",
             "sstart", "send", "evalue", "bitscore"]
    # Extract scores
    df2 = pd.read_csv(child, sep="\t", names=names, index_col=0)
    s = df2.iloc[:, 10]
    # Set Series name
    s.name = proteome
    # Add scores to the df
    if len(s) > 1000:
        df = pd.concat([df, s], axis=1)

# Set index name
df.index.names = ["Genes"]
# Fill missing values
df = df.fillna(1)
# Lenght normalization
gene_norm = df["UP000000803_7227"]
df_gnorm = df.divide(gene_norm, axis="index")
# Species normalization
df_norm = df_gnorm.subtract(df_gnorm.mean(axis=1), axis=0).divide(df_gnorm.std(axis=1), axis=0)
#for column in df.columns:
#    df[column] = (df[column] - df[column].mean()) / df[column].std()

print(df_norm)
# Save df
df_norm.to_csv("temporal")
