#!/usr/bin/env python3

"Converts a fasta nucleotide sequence into its corresponding protein translation"

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click


def translate_cds(seqs, output):
    "Converts a fasta nucleotide sequence into its corresponding protein translation"
    
    # Load coding sequences
    seqs = SeqIO.parse(seqs, 'fasta')
    # Translate CDS and store them into a list
    prots = [SeqRecord(seq.seq.translate(to_stop=True), seq.id, seq.name, seq.description) for seq in seqs]
    # Save to file if a name is provided
    if output is not None:
        SeqIO.write(prots, output, 'fasta')
    # Print to standard output if no name is given
    if output is None:
        for prot in prots:
            # Fasta format style
            print(f'>{prot.description}\n{prot.seq}')

# Command line interface (CLI) options
SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=SETTINGS)
@click.option('-s',
              '--seqs',
              help='Coding sequences to translate')
@click.option('-o',
              '--output',
              help='Output file name. If nothing is provided, output is sent to STDOUT')
 
# CLI function
def cli(seqs, output):
    translate_cds(seqs, output)

if __name__ == '__main__':
    cli()
