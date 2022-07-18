rule all:
    input:
        "diptera_orthogroups.csv"

rule extract_orthogroups:
    input:
        ogs="../resources/odb10v1_OGs.tab",
        og2genes="../resources/odb10v1_OG2genes.tab",
        genes="../resources/odb10v1_genes.tab",
        species="../resources/odb10v1_species.tab",
        mappings="../resources/orthodb2flybase_table.csv"

    output:
        "diptera_orthogroups.csv"
    shell:
        """
        scripts/filter_orthodb_tables.py --orthogroups {input.ogs} --taxid 7147 --og2genes {input.og2genes} --genes {input.genes} --species {input.species} -m {input.mappings} -o {output}
        """
