#!/usr/bin/env python3

"""
Returns the longest isoform of each CDS based on its gene=[ID] feature
"""

from Bio import SeqIO
import click


def get_longest_isoform(seqs: SeqIO.SeqRecord) -> str:
    """
    Returns longest isoform of each protein from a multifasta
    """

    # Read initial sequences into a SeqIO object
    seqs = SeqIO.parse(seqs, 'fasta')
    # Initialize a dict to store filtered sequences
    fseqs = {}
    # Iterate over all sequences
    for seq in seqs:
        # Get gene feature from NCBI CDS fasta header
        seq_id = seq.description.split()[1]
        # If the gene is absent on the filtered dict, add it
        if seq_id not in fseqs.keys():
            fseqs[seq_id] = seq
        # If the gene is already present, replace it with the new sequence if it is longer than the existing one
        elif len(seq) > len(fseqs[seq_id]):
            fseqs[seq_id] = seq

    return fseqs


def save_cds_sequences(sequences, organism):

    # Get longest isoforms
    fseqs = get_longest_isoform(sequences)
    # Set output name
    output = organism + "_filt.fna"
    # Save filtered sequences into output
    with open(output, "w") as fh:
        for key in fseqs:
            seq = fseqs[key]
            fh.write(f'>{seq.id}\n{seq.seq}\n')


def save_protein_sequences(sequences, organism):

    # Get longest isoforms
    fseqs = get_longest_isoform(sequences)
    # Set output name
    output = organism + "_filt.faa"
    # Save filtered sequences into output
    with open(output, "w") as fh:
        for key in fseqs:
            seq = fseqs[key]
            # Translate CDS into amino acid sequence
            seq.seq = seq.seq.translate(to_stop=True)
            fh.write(f'>{seq.id}\n{seq.seq}\n')


# Set CLI parameters
SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=SETTINGS)
@click.option("-s",
              "--sequences",
              help="Sequences to filter")
@click.option("-o",
              "--organism",
              help="Organism name used for output")
# Command line interface
def cli(sequences, organism):
    """
    Returns the longest isoform of each CDS based on its gene=[ID] feature
    """

    save_cds_sequences(sequences, organism)
    # save_protein_sequences(sequences, organism)


if __name__ == "__main__":
    cli()
