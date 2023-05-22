import pandas as pd
SPECIES = pd.read_csv("../resources/drosophila_urls.csv", names=["organism", "url"]).organism
ORTHOGROUPS = pd.read_csv("../resources/Orthogroups.tsv", sep="\t", usecols=["Orthogroup"]).Orthogroup

rule all:
    input:
        expand("../results/{species}_filt.fna", species=SPECIES),
        expand("../results/evolutionary_rates/{og}.aln.fna.FUBAR.json", og=ORTHOGROUPS),
        expand("../results/evolutionary_rates/{og}.aln.fna.ABSREL.json", og=ORTHOGROUPS)


rule download_genomes:
    input:  "../resources/drosophila_urls.csv"
    output: temp(expand("../resources/genomes/{species}.fa.gz", species=SPECIES))
    shell:
        """
        scripts/download_ncbi_genomes.py --list {input} --folder ../resources/genomes/
        """


rule decompress_genomes:
    input:  "../resources/genomes/{species}.fa.gz"
    output: temp("../results/{species}.fa")
    shell:
        """
        gunzip -c {input} > {output}
        """


rule get_longest_isoform:
    input:  "../results/{species}.fa"
    output: temp("../results/{species}_filt.fna")
    shell:
        """
        ./scripts/get_longest_isoform.py -s {input} -o {wildcards.species} && \
        mv {wildcards.species}_filt.fna ../results/
        """


rule extract_orthogroup_multifasta:
    input: expand("../results/{species}_filt.fna", species=SPECIES)
    params: "../results/evolutionary_rates/"
    output: expand("../results/evolutionary_rates/{og}.fna", og=ORTHOGROUPS),
            expand("../results/evolutionary_rates/{og}.faa", og=ORTHOGROUPS)
    shell:
        """
        cat {input} > sequences.fna && \
        scripts/extract_orthogroups.py -o ../resources/Orthogroups.tsv -s sequences.fna -d {params} && \
        rm sequences.fna
        """


rule protein_alignment:
    input: "../results/evolutionary_rates/{og}.faa"
    output: "../results/evolutionary_rates/{og}.aln.faa"
    shell:
      """
      mafft {input} > {output}
      """

rule codon_alignment:
    input: "../results/evolutionary_rates/{og}.aln.faa",
           "../results/evolutionary_rates/{og}.fna"
    output: "../results/evolutionary_rates/{og}.aln.fna"
    shell:
      """
      pal2nal.pl {input} -output fasta > {output}
      """


rule run_fubar_analysis:
    input: "../results/evolutionary_rates/{og}.aln.fna",
           "../resources/SpeciesTree_rooted.txt"
    output: "../results/evolutionary_rates/{og}.aln.fna.FUBAR.json"
    shell:
      """
      hyphy fubar --alignment {input[0]} --tree {input[1]}
      """


rule run_absrel_analysis:
    input: "../results/evolutionary_rates/{og}.aln.fna",
           "../resources/SpeciesTree_rooted.txt"
    output: "../results/evolutionary_rates/{og}.aln.fna.ABSREL.json"
    shell:
      """
      hyphy absrel --alignment {input[0]} --tree {input[1]}
      """
