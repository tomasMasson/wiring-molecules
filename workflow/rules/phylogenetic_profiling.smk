from pathlin import Path
path = Path("../resources/uniprot_proteomes/")
PROTS = [p.stem.split(".")[0] for p in path.iterdir()]


rule all:
    input:
        "../results/phylogenetic_profiles/conservation_matrix.csv"


rule decompress_proteomes:
    input:
        "../resources/uniprot_proteomes/{proteomes}.fasta.gz"
    output:
        temp("../resources/uniprot_proteomes/{proteomes}.fasta")
    shell:
        """
        gunzip --keep {input}
        """


rule add_species_names:
    input:
        "../resources/uniprot_proteomes/{proteome}.fasta",
    output:
        temp("../resources/uniprot_proteomes/{proteome}_renamed.fasta")
    shell:
        """
        ./scripts/add_species_names.py --input {input} > {output}
        """


rule make_diamond_db:
    input:
        "../resources/uniprot_proteomes/{proteome}_renamed.fasta"
    output:
        "../results/phylogenetic_profiles/{proteome}.dmnd"
    shell:
        "diamond makedb --in {input} --db {output}"


rule diamond_conservation_scoring:
    input:
        proteome="../resources/uniprot_proteomes/UP000000803_7227_renamed.fasta",
        db="../results/phylogenetic_profiles/{proteome}.dmnd"
    output:
        "../results/phylogenetic_profiles/{proteome}_conservation.tsv"
    shell:
        "diamond blastp --query {input.proteome} --db {input.db} --max-target-seqs 1 --ultra-sensitive --outfmt 6 --out {output}"


rule aggregate_diamond_scoring:
    input:
        expand("../results/phylogenetic_profiles/{proteome}_conservation.tsv", proteome=PROTS)
    output:
        "../results/phylogenetic_profiles/diamond_searches.tsv"
    shell:
        "cat {input} >> {output}"


rule build_conservation_matrix:
    input:
        "../results/phylogenetic_profiles/diamond_searches.tsv"
    output:
        "../results/phylogenetic_profiles/conservation_matrix.csv"
    shell:
        "./scripts/build_conservation_matrix.py --input {input} --output {output}"
