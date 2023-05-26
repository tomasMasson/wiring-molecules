import pandas as pd

SPECIES = pd.read_csv("../resources/drosophila_urls.csv", names=["organism", "url"]).organism


rule:
    input:
        "../results/evo_rates/Orthogroups.tsv",
        "../results/evo_rates/SpeciesTree.txt",
        "../results/evo_rates/cds_database.fna"


rule download_cds_data:
    input:  "../resources/drosophila_urls.csv"
    output: temp(expand("{species}.fa.gz", species=SPECIES))
    shell: "scripts/download_ncbi_genomes.py --list {input}"


rule decompress_cds_data:
    input:  "{species}.fa.gz"
    output: temp("{species}.fa")
    shell: "gunzip -c {input} > {output}"


rule get_longest_isoform:
    input:  "{species}.fa"
    output: temp("{species}.fna")
    shell: "./scripts/get_longest_isoform.py -s {input} -o {output}"


rule translate_proteomes:
    input: "{species}.fna"
    output: temp("{species}.faa")
    shell: "scripts/translate_cds.py -s {input} > {output}"


rule create_orthofinder_input:
    input: expand("{species}.faa", species=SPECIES)
    output: directory("proteome/")
    shell: "mkdir {output} && mv {input} {output}"


rule run_orthofinder:
    input: "proteome"
    params: outdir="of",
            name="of"
    output:
        "../results/evo_rates/Orthogroups.tsv",
        "../results/evo_rates/SpeciesTree.txt",
    shell:
        """
        orthofinder -f {input} -o {params.outdir} -n {params.name} && \
        cp {params.outdir}/Results_{params.name}/Orthogroups/Orthogroups.tsv {output[0]} && \
        cp {params.outdir}/Results_{params.name}/Species_Tree/SpeciesTree_rooted.txt {output[1]} && \
        rm {input} {params.outdir} -rf 
        """


rule build_cds_database:
    input:  expand("{species}.fna", species=SPECIES)
    output: "../results/evo_rates/cds_database.fna"
    shell: "cat {input} > {output}"
