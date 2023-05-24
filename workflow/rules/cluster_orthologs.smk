import pandas as pd

SPECIES = pd.read_csv("../resources/drosophila_urls.csv", names=["organism", "url"]).organism


rule:
    input:
        "../results/orthofinder/Orthogroups.tsv"


rule download_proteome_data:
    input:  "../resources/drosophila_urls.csv"
    output: temp(expand("{species}.fa.gz", species=SPECIES))
    shell: "scripts/download_ncbi_genomes.py --list {input}"


rule decompress_proteome_data:
    input:  "{species}.fa.gz"
    output: temp("{species}.fa")
    shell: "gunzip -c {input} > {output}"


rule get_longest_proteoform:
    input:  "{species}.fa"
    output: "{species}.fasta"
    shell: "./scripts/get_longest_isoform.py -s {input} -o {wildcards.species}"


rule create_orthofinder_input:
    input: expand("{species}.fasta", species=SPECIES)
    output: directory("proteome/")
    shell: "mkdir {output} && mv {input} {output}"


rule run_orthofinder:
    input: "proteome"
    params: outdir="of",
            name="of"
    output: "../results/orthofinder/Orthogroups.tsv"
    shell:
        """
        orthofinder -f {input} -o {params.outdir} -n {params.name} && \
        cp {params.outdir}/Results_{params.name}/Orthogroups/Orthogroups.tsv {output}
        rm {params.outdir} -rf
        """
