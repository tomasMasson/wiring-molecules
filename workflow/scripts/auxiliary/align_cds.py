#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import click
import subprocess

def translate_cds(sequences, output):
    "Converts a fasta nucleotide sequence into its corresponding protein translation"

    # Load coding sequences
    seqs = SeqIO.parse(sequences, 'fasta')
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

def align_proteins(proteins, output):
    subprocess.run(f'mafft {proteins} > {output}', shell=True)

def align_cds(prot_aln, seqs, output):
    subprocess.run(f'pal2nal.pl {prot_aln} {seqs} -output paml > {output}', shell=True)

@click.command()
@click.option('-d',
              'directory',
              help='Directory with coding sequences')
def cli(directory):
    p = Path(directory)
    for file in p.iterdir():
        prot_file = file.parent / f'{file.stem}.faa'
        prot_aln = file.parent / f'{file.stem}.aln.faa'
        cds_aln = file.parent / f'{file.stem}.aln.fna'
        translate_cds(file, prot_file)
        align_proteins(prot_file, prot_aln)
        align_cds(prot_aln, file, cds_aln)



if __name__ == '__main__':
    cli()
