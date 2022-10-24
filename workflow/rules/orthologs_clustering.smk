from pathlib import Path

GENOMES, = glob_wildcards("../resources/{genome}.all.fa")

rule all:
    input: 
        "../results/orthology_inference/Log.txt",
        "../results/orthology_inference/proteome_cds_database.fna"

rule download_flyxcdb_data:
    params:
        "http://prodata.swmed.edu/FlyXCDB/info.list.new21_26.html"
    output:
        "../resources/flyxcdb_data.csv"
    conda:
        "envs/wiring_evolution.yml"
    shell:
        """
        scripts/download_flyxcdb_data.py --html {params} --output {output}
        """

rule filter_longest_gene_isoform:
    input:
        "../resources/{genome}.all.fa"
    output:
        temp("../results/orthology_inference/{genome}.all.fa")
    shell:
        "scripts/filter_longest_isoform.py --seqs {input} > {output}"

rule aggregate_coding_sequences:
    input:
        expand("../results/orthology_inference/{genome}.all.fa", genome=GENOMES)
    output:
        "../results/orthology_inference/proteome_cds_database.fna"
    shell:
        "cat {input} > {output}"

rule translate_cds_proteomes:
    input:
        "../results/orthology_inference/{genome}.all.fa"
    output:
        "../results/orthology_inference/{genome}.faa"
    shell:
        "scripts/translate_cds.py --seqs {input} --output {output}"

rule run_orthofinder:
    input:
        expand("../results/orthology_inference/{genome}.faa", genome=GENOMES)
    params:
        "../results/orthology_inference/"
    output:
        "../results/orthology_inference/Log.txt"
    shell:
        """
        orthofinder -f {params} -n of > {output} && \
        mv {params}/OrthoFinder/Results_of/* {params} && \
        rm -rf {params}OrthoFinder/
        """
