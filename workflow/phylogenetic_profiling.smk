PROTS, = glob_wildcards("../resources/{proteome}.fasta")

rule all:
    input:
        "../results/phylogenetic_profiles/diamond_results.tsv"

rule make_diamond_db:
    input:
        "../resources/UP000000803_7227.fasta"
    output:
        "../results/phylogenetic_profiles/diamond_db.dmnd"
    shell:
        "diamond makedb --in {input} --db {output}"

rule diamond_scoring:
    input:
        proteome="../resources/{PROTS}.fasta",
        db="../results/phylogenetic_profiles/diamond_db.dmnd"
    output:
        "../results/phylogenetic_profiles/{PROTS}_results.tsv"
    shell:
        "diamond blastp --query {input.proteome} --db {input.db} --max-target-seqs 1 --ultra-sensitive --outfmt 6 --out {output}"

rule aggregate_diamond_scoring:
    input:
        expand("../results/phylogenetic_profiles/{proteomes}_results.tsv", proteomes=PROTS)
    output:
        "{OUTDIR}diamond_results.tsv"
    shell:
        "cat {input} >> {output}"
