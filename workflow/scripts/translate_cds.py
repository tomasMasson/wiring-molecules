#!/usr/bin/env python3

"Converts a fasta nucleotide sequence into its corresponding protein translation"

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click


def translate_cds(seq: SeqIO.SeqRecord) -> SeqIO.SeqRecord:
    "Returns a CDS translation in the standard or invertebrate mitochondrial codes"
    # Get protein sequence (NCBI translation table=1)
    t = seq.seq.translate(to_stop=True)
    # Check that CDS and protein lengths matches
    if len(t) == len(seq.seq)/3 - 1:
        return SeqRecord(t, seq.id, seq.name, seq.description)
    # Treat the CDS as mitochondrial (table=5)
    else:
        t = seq.seq.translate(to_stop=True, table=5)
        return SeqRecord(t, seq.id, seq.name, seq.description)


def translate_multifasta(seqs: SeqIO.SeqRecord, output: str) -> None:
    "Converts a CDS multifasta into its protein version"
    # Load coding sequences
    sequences = SeqIO.parse(seqs, 'fasta')
    # Translate CDS and store them into a list
    prots = [translate_cds(seq)
             for seq in sequences
             if len(seq.seq) % 3 == 0]
    # Save to file if a name is provided
    if output is not None:
        SeqIO.write(prots, output, 'fasta')
    # Print to standard output if no name is given
    if output is None:
        for prot in prots:
            # Fasta format style
            print(f'>{prot.description}\n{prot.seq}')
    return None

# Command line interface (CLI) options
SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=SETTINGS)
@click.option('-s',
              '--seqs',
              help='Coding sequences to translate')
@click.option('-o',
              '--output',
              help='Output file name. If not provided, output is sent to STDOUT')
 
# CLI function
def cli(seqs, output):
    translate_multifasta(seqs, output)

if __name__ == '__main__':
    cli()
