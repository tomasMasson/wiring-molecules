import pandas as pd


def filter_scorthogroups(data: str, outfile=None) -> pd.Series:
    "Keep only single copy Orthogroups present in all species"
    # Remove non-universal orthogroups
    df = pd.read_csv(data, sep="\t").dropna()
    # Collect Orthogroup names that are not single copy
    exclude = []
    for row in df.iterrows():
        for i in row[1]:
            if len(i.split(",")) > 1:
                exclude.append(row[1].Orthogroup)  
                break
    # Save scOrthogroups
    dff = df[~df.Orthogroup.isin(exclude)]
    if outfile:
        dff.to_csv(outfile, index=False)
    return dff.Orthogroup


SPECIES = pd.read_csv("../resources/drosophila_urls.csv", names=["organism", "url"]).organism
ORTHOGROUPS = filter_scorthogroups("../results/evo_rates/Orthogroups.tsv")


rule all:
    input:
        expand("../results/evo_rates/{og}.aln.fna.FUBAR.json", og=ORTHOGROUPS),
        expand("../results/evo_rates/{og}.aln.fna.ABSREL.json", og=ORTHOGROUPS)


rule filter_single_copy_ortogroups:
    input: "../results/evo_rates/Orthogroups.tsv"
    output: temp("../results/evo_rates/Orthogroups_filtered.tsv")
    run:
        filter_scorthogroups(input[0], output[0])


rule extract_orthogroup_multifasta:
    input:
        "../results/evo_rates/Orthogroups_filtered.tsv",
        "../results/evo_rates/cds_database.fna"
    params: "../results/evo_rates/"
    output: expand("../results/evo_rates/{og}.fna", og=ORTHOGROUPS),
            expand("../results/evo_rates/{og}.faa", og=ORTHOGROUPS)
    shell:
        "scripts/extract_orthogroups.py -o {input[0]} -s {input[1]} -d {params}"


rule protein_alignment:
    input: "../results/evo_rates/{og}.faa"
    output: "../results/evo_rates/{og}.aln.faa"
    shell:
      """
      mafft {input} > {output}
      """

rule codon_alignment:
    input: "../results/evo_rates/{og}.aln.faa",
           "../results/evo_rates/{og}.fna"
    output: "../results/evo_rates/{og}.aln.fna"
    shell:
      """
      pal2nal.pl {input} -output fasta > {output}
      """


rule run_fubar_analysis:
    input: "../results/evo_rates/{og}.aln.fna",
           "../results/evo_rates/SpeciesTree.txt"
    output: "../results/evo_rates/{og}.aln.fna.FUBAR.json"
    shell:
      """
      hyphy fubar --alignment {input[0]} --tree {input[1]}
      """


rule run_absrel_analysis:
    input: "../results/evo_rates/{og}.aln.fna",
           "../results/evo_rates/SpeciesTree.txt"
    output: "../results/evo_rates/{og}.aln.fna.ABSREL.json"
    shell:
      """
      hyphy absrel --alignment {input[0]} --tree {input[1]}
      """
