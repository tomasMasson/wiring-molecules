PROTS, = glob_wildcards("../results/phylogenetic_profiles/test_proteomes/{proteome}.fasta")

rule all:
    input:
        "../results/phylogenetic_profiles/conservation_matrix.test.csv"


rule add_species_names:
    input:
        "../results/phylogenetic_profiles/test_proteomes/{proteome}.fasta",
    output:
        "../results/phylogenetic_profiles/test_proteomes/{proteome}_renamed.fasta"
    shell:
        "./scripts/add_species_names.py --input {input} > {output} "

rule make_diamond_db:
    input:
        "../results/phylogenetic_profiles/test_proteomes/{proteome}_renamed.fasta"
    output:
        "../results/phylogenetic_profiles/{proteome}.test.dmnd"
    shell:
        "diamond makedb --in {input} --db {output}"


rule diamond_conservation_scoring:
    input:
        proteome="../results/phylogenetic_profiles/test_proteomes/UP000000803_7227_renamed.fasta",
        db="../results/phylogenetic_profiles/{proteome}.test.dmnd"
    output:
        "../results/phylogenetic_profiles/{proteome}_conservation.test.tsv"
    shell:
        "diamond blastp --query {input.proteome} --db {input.db} --max-target-seqs 1 --ultra-sensitive --outfmt 6 --out {output}"


rule aggregate_diamond_scoring:
    input:
        expand("../results/phylogenetic_profiles/{proteome}_conservation.test.tsv", proteome=PROTS)
    output:
        "../results/phylogenetic_profiles/diamond_searches.test.tsv"
    shell:
        "cat {input} >> {output}"


rule build_conservation_matrix:
    input:
        "../results/phylogenetic_profiles/diamond_searches.test.tsv"
    output:
        "../results/phylogenetic_profiles/conservation_matrix.test.csv"
    shell:
        "./scripts/build_conservation_matrix.py --input {input} --output {output}"
