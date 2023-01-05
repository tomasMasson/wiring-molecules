SAMPLES = ["t4_1_hrp_1",
           "t4_1_hrp_2",
           "t4_2_hrp_1",
           "t4_2_hrp_2",
           "t4_1_h2o2_1",
           "t4_1_h2o2_2",
           "t4_2_h2o2_1",
           "t4_2_h2o2_2",
           "t5_1_hrp_1",
           "t5_1_hrp_2",
           "t5_2_hrp_1",
           "t5_2_hrp_2",
           "t5_1_h2o2_1",
           "t5_1_h2o2_2",
           "t5_2_h2o2_1",
           "t5_2_h2o2_2"]
FORMATS = ["csv", "png"]

rule all:
    input: 
        "t4_consensus_surfaceome.csv", 
        "t5_consensus_surfaceome.csv",
        "pca_analysis.png",
        "venn_diagram.png"

rule run_analysis:
    input:
        dataset="../resources/t4vt5_pl_dataset.csv",
        mappings="../resources/uniprot2flybase.tab",
        go_terms="../resources/gene_association.fb",
        flyxcdb="../resources/flyxcdb_data.csv"
    params:
        sample="{sample}"
    output:
        "{sample}.{ext}"
    shell:
        "scripts/ratiometric_classifier.py --data {input.dataset} --mappings {input.mappings} --go_terms {input.go_terms} --flyxcdb {input.flyxcdb} --label {params.sample}"

rule merge_outputs:
    input:
        expand("{sample}.csv", sample=SAMPLES)
    output:
        temp("t4vt5_signficant_proteins.csv")
    shell:
        """
        cat {input} > {output} && \
        rm t4_*.csv t5_*.csv
        """

rule get_consensus_surfaceomes:
    input:
        "t4vt5_signficant_proteins.csv"
    output:
        "t4_consensus_surfaceome.csv", 
        "t5_consensus_surfaceome.csv"
    shell:
        """
        scripts/get_consensus_surfaceomes.py --input {input}
        """

rule plot_pca_analysis:
    input:
        "../resources/t4vt5_pl_dataset.csv"
    output:
        "pca_analysis.png"
    shell:
        "scripts/plot_pca.py --data {input}"

rule plot_venn_diagram:
    input:
        "t4_consensus_surfaceome.csv",
        "t5_consensus_surfaceome.csv"
    output:
        "venn_diagram.png"
    shell:
        "scripts/plot_venn_diagram.py --data1 {input[0]} --data2 {input[1]}"
