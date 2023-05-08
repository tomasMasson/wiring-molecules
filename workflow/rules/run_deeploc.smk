from Bio import SeqIO
from pathlib import Path
import textwrap
import subprocess


RANGE = range(0, config["PARAMS"]["size"], config["PARAMS"]["chunk"])


def write_seq_object(seq: SeqIO.SeqRecord) -> str:
    # Cap line length to 60 characters
    sequence = '\n'.join(textwrap.wrap(str(seq.seq), 60))
    # Save output 
    return f">{seq.description}\n{sequence}\n"


def split_multifasta(proteome: str) -> None:
    "Divides a multifasta file into smaller ones based on the chuck value"
    # Load sequences
    seqs = list(SeqIO.parse(proteome, "fasta"))
    # Split multifasta into chunks of 1000
    for i in range(0, len(seqs), 1000):
      with open(f"proteome_{i}.faa", "w") as fh:
        if i < len(seqs):
          for seq in seqs[i: i+1000]:
              fh.write(write_seq_object(seq))
        else:
          for seq in seqs[i:]:
              fh.write(write_seq_object(seq))
    return None


rule all:
    input: 
        "deeploc2_predictions.csv"


rule download_proteome:
    params:
        config["DATA"]["proteome_url"]
    output:
        config["DATA"]["proteome"]
        # temp(config["DATA"]["proteome"])
    shell:
        """
        wget {params} -O tmp.gz && \
        gunzip -c tmp.gz > {output} && \
        rm tmp.gz 
        """


rule split_multifasta:
    input:
        config["DATA"]["proteome"]
    output:
        expand("proteome_{chunk}.faa", chunk=RANGE)
    run:
        split_multifasta(input[0])


rule run_deeploc2:
    input: 
        expand("proteome_{chunk}.faa", chunk=RANGE)
    params:
        "deeploc2_output"
    output:
        "deeploc2_predictions.csv"
    run:
        for i in input:
            subprocess.run(["deeploc2", "-f", i, "-o", params[0]])
        with open(output[0], "w") as fh:
            p = Path(params[0])
            for file in p.iterdir():
                with open(file) as infile:
                    fh.write(infile.read())


