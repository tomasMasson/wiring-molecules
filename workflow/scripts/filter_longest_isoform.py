#!/usr/bin/env python3

"""
Retain the longest isoform of each coding sequence, either nucleotide or protein, based on the gene:[ID] attribute from the fasta header.
"""

from typing import NewType
from Bio import SeqIO
import click


def get_longest_isoform(seqs: SeqIO.SeqRecord) -> str:
    """
    Print the longest isoform for each protein to standard output 
    """

    # Read initial sequences into a SeqIO object
    seqs = SeqIO.parse(seqs, 'fasta')
    # Initialize a dict to store filtered sequences
    seqs_filt = {}
    # Iterate over all sequences
    for seq in seqs:
        # Get the gene feature on fasta header 
        seq_id = seq.description.split()[3]
        # If the header is not present in the filtered dict, just add the sequence
        if seq_id not in seqs_filt.keys():
            seqs_filt[seq_id] = seq
        # If the gene is already present, add the new sequence if it is longer than the existing one 
        elif len(seq) > len(seqs_filt[seq_id]):
            seqs_filt[seq_id] = seq
    # Send filtered sequences into standard output
    for key in seqs_filt:
        seq = seqs_filt[key]
        print(f'>{seq.description}\n{seq.seq}')

# Set CLI parameters
## Add help page if no variable is provided through th command line
SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=SETTINGS)
@click.option("-s",
              "--seqs",
              "sequences",
              help="Sequences to filter")

# Put all together into a principal function
def cli(sequences):
    "Command line interface"
    get_longest_isoform(sequences)


if __name__ == "__main__":
    cli()
