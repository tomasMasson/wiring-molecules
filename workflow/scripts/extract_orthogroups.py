#!/usr/bin/env python3.8

"""
"""


import click
import pandas as pd
from Bio import SeqIO


def extract_orthogroups(orthogroups, multifasta, outdir):

    seqs = SeqIO.to_dict(SeqIO.parse(multifasta, "fasta"))
    df = pd.read_csv(orthogroups, sep="\t").dropna()
    dc = {}
    for row in df.iterrows():
        og = row[1].values[0]
        species = row[1].index[1:]
        genes = row[1].values[1:]
        orthologs = [(specie, gene)
                     for specie, gene in zip(species, genes)
                     ]
        p = outdir + f"{og}.fna"
        with open(p, "w") as fh:
            for seq in orthologs:
                seq_id = seq[0]
                seq_seq = seqs[seq[1]].seq
                fh.write(f">{seq_id}\n{seq_seq}\n")
        p = outdir + f"{og}.faa"
        with open(p, "w") as fh:
            for seq in orthologs:
                seq_id = seq[0]
                seq_seq = seqs[seq[1]].seq
                fh.write(f">{seq_id}\n{seq_seq.translate(to_stop=True)}\n")
        # p = outdir + f"{key}.faa"
        # with open(p, "w") as fh:
        #     for i in dc[key]:
        #         seq_data = seqs[i]
        #         seq_id = seq_data.id.split("_cds")[0].split("|")[1]
        #         seq_seq = seq_data.seq
        #         fh.write(f">{seq_id}\n{seq_seq.translate(to_stop=True)}\n")
        # key = row[1].values[0]
        # genes = []
        # for column in row[1].values[1:]:
        #     for value in column.split(", "):
        #         genes.append(value)
        # if len(genes) == 6:
        #     dc[key] = genes
    # for key in dc:
        # p = outdir + f"{key}.fna"
        # with open(p, "w") as fh:
        #     for i in dc[key]:
        #         seq_data = seqs[i]
        #         seq_id = seq_data.id.split("_cds")[0].split("|")[1]
        #         seq_seq = seq_data.seq
        #         fh.write(f">{seq_id}\n{seq_seq}\n")
        # p = outdir + f"{key}.faa"
        # with open(p, "w") as fh:
        #     for i in dc[key]:
        #         seq_data = seqs[i]
        #         seq_id = seq_data.id.split("_cds")[0].split("|")[1]
        #         seq_seq = seq_data.seq
        #         fh.write(f">{seq_id}\n{seq_seq.translate(to_stop=True)}\n")


# CLI options
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command()
@click.option("-o",
              "--orthogroups",
              help="")
@click.option("-s",
              "--sequences",
              help="")
@click.option("-d",
              "--directory",
              help="")
def cli(orthogroups, sequences, directory):
    """
    """
    extract_orthogroups(orthogroups, sequences, directory)


if __name__ == "__main__":
    cli()
