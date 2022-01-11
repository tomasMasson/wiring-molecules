#!/usr/bin/env python3

"""
Keeps the longest isoform of each coding sequence (CDS) in multifasta downloaded from FlyBase genome FTP.
"""

from typing import NewType
from Bio import SeqIO
import click


def print_sequence_dict(seqs: dict) -> str:
    """
    Prints a multi sequence dictionary in fasta format
    """

    # Iterate over the sequences dictionary
    for seq in seqs:
        sequence = seqs[seq]
        # Print sequence in FASTA format
        print(f">{sequence.description}\n{sequence.seq}")


def get_identifier(seq: SeqIO.SeqRecord) -> str:
    """
    Parse the Flybase gene identifier from the fasta header. For this, it uses it the description attribute of SeqRecord object from Biopython SeqIO module.
    In this particular case, the parser is prepared to deal with the fasta header implemented on genomic CDS multifasta files (e.g. http://ftp.flybase.net/genomes/Drosophila_willistoni/dwil_r1.3_FB2015_01/fasta).
    """

    gene_name = seq.description.split(';')[6][8:19]
    return gene_name

def compare_length(seq: SeqIO.SeqRecord) -> str:
    """
    """

    length = len(seq) 

    return length

def get_longest_isoform(seqs: SeqIO.SeqRecord) -> str:
    """
    Print to standard output the longest isoform for each protein.
    """

    # Read initial sequences into a SeqIO object
    seqs = SeqIO.parse(seqs, 'fasta')
    # Initialize a sequence Dict, using faste header as key
    seqs_dict = SeqIO.to_dict(seqs,
                              key_function=lambda rec: rec.description)
    # Initialize a dict to store filtered sequences
    seqs_filt = {}
    # For each sequence in the dict
    for key in seqs_dict:
        # Retrieve the Seq object
        seq = seqs_dict[key]
        # Get sequence identifier
        seq_id = get_identifier(seq)
        # If the identifier is new in the filtered dict, just add the sequence
        if seq_id not in seqs_filt.keys():
            seqs_filt[seq_id] = seq
        # The other option is that the identifier is already on the filtered dict.
        # In that case, add the new sequence if it is longer than the existing one 
        elif len(seq) > len(seqs_filt[seq_id]):
            seqs_filt[seq_id] = seq
    # Sent the filtered sequence dict into standard output
    print_sequence_dict(seqs_filt)


# Set CLI parameters
## Add help page if no variable is provided through th command line
SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=SETTINGS)
@click.option("-s",
              "--seqs",
              "sequences",
              help="Sequences to filtered")

# Put all together into a principal function
def cli(sequences):
    get_longest_isoform(sequences)


if __name__ == "__main__":
    cli()
