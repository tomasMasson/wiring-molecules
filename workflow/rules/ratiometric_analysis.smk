import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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

# FORMATS = ["csv", "png"]


rule all:
    input:
        "../results/t4t5_surface_proteomics/t4t5_consensus_surfaceome.csv",
        "../results/t4t5_surface_proteomics/pca_analysis.svg",
        "../results/t4t5_surface_proteomics/venn_diagram.svg",
        "../results/t4t5_surface_proteomics/t4t5_surfaceome_physical_interactions.csv",
        "../results/t4t5_surface_proteomics/cell_adhesion_molecules_abundance.svg"


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
        """
        scripts/ratiometric_classifier.py --data {input.dataset} --mappings {input.mappings} --go_terms {input.go_terms} --flyxcdb {input.flyxcdb} --label {params.sample}
        """


rule merge_outputs:
    input:
        expand("{sample}.csv", sample=SAMPLES)
    output:
        temp("t4t5_signficant_proteins.csv")
    shell:
        # cat {input} > {output} && \
        # rm t4_*.csv t5_*.csv
        """
        cat {input} > {output}
        """


rule get_consensus_surfaceomes:
    input:
        "t4t5_signficant_proteins.csv"
    output:
        "t4t5_consensus_surfaceome.csv"
    shell:
        """
        scripts/get_consensus_surfaceomes.py --input {input}
        """


rule plot_pca_analysis:
    input:
        "../resources/t4vt5_pl_dataset.csv"
    output:
        "pca_analysis.svg"
    shell:
        """
        scripts/plot_pca.py --data {input}
        """


rule plot_venn_diagram:
    input:
        "t4t5_consensus_surfaceome.csv"
    output:
        "venn_diagram.svg"
    shell:
        """
        scripts/plot_venn_diagram.py --data {input}
        """


rule move_outputs:
    input:
        "t4t5_consensus_surfaceome.csv",
        "pca_analysis.svg",
        "venn_diagram.svg"
    params:
        "../results/t4t5_surface_proteomics/"
    output:
        "../results/t4t5_surface_proteomics/t4t5_consensus_surfaceome.csv",
        "../results/t4t5_surface_proteomics/pca_analysis.svg",
        "../results/t4t5_surface_proteomics/venn_diagram.svg"
    shell:
        """
        mv {input} t4_* t5_* {params}
        """

# rule get_flybase_physical_interactions:
#     params: "https://ftp.flybase.net/releases/current/precomputed_files/genes/physical_interactions_mitab_fb_2024_06.tsv.gz"
#     output:
#         "../results/t4t5_surface_proteomics/flybase_physical_interactions.tsv.gz"
#     shell:
#         "wget {params} -O {output}"

rule get_surfaceome_physical_interactions:
    input:
        "../results/t4t5_surface_proteomics/t4t5_consensus_surfaceome.csv",
        "../resources/physical_interactions_mitab_fb_2025_03.tsv.gz",
        "../resources/flyxcdb_data.csv"
    output:
        "../results/t4t5_surface_proteomics/t4t5_surfaceome_physical_interactions.csv"
    run:
        df = pd.read_csv(input[0])
        interactions = pd.read_csv(input[1], sep="\t")
        flyxcdb = pd.read_csv(input[2])
        interactions["#ID(s) Interactor A"] = interactions["#ID(s) Interactor A"].str.split(":", expand=True)[1]
        interactions["ID(s) Interactor B"] = interactions["ID(s) Interactor B"].str.split(":", expand=True)[1]
        subset = interactions[interactions["#ID(s) Interactor A"].isin(df.Protein)]
        subset2 = subset[subset["ID(s) Interactor B"].isin(df.Protein)]
        subset3 = subset2[subset2["#ID(s) Interactor A"].isin(flyxcdb.GeneID)]
        subset4 = subset3[subset3["ID(s) Interactor B"].isin(flyxcdb.GeneID)]
        subset4.to_csv(output[0])

rule plot_wiring_proteins_abundances:
    input:
        "../results/t4t5_surface_proteomics/t4t5_consensus_surfaceome.csv",
        "../resources/expanded_cams_list.csv"
    output:
        "../results/t4t5_surface_proteomics/cell_adhesion_molecules_abundance.svg",
    run:
        # Load proteomic data
        df = pd.read_csv(input[0])
        # Load cell adhesion molecules data
        cam = pd.read_csv(input[1])
        # Annotate proteins as CAMs
        df["CAM"] = df["Protein"].isin(cam["FBgn"])
        sns.set(font_scale=1.2)
        sns.set_style("white")
        fig, ax = plt.subplots(figsize=(5,5))
        dff = df.replace({True: "T4 wiring",
                          False: "T4 surfaceome"})
        sns.boxplot(x="CAM",
                    y="T4_signal",
                    data=dff[dff["T4_signal"] != 0],
                    color="#f1a340ff",
                    showcaps=False,
                    linewidth=2,
                    ax=ax)
        dff = df.replace({True: "T5 wiring",
                          False: "T5 surfaceome"})
        sns.boxplot(x="CAM",
                    y="T5_signal",
                    data=dff[dff["T5_signal"] != 0],
                    color="#998ec3ff",
                    showcaps=False,
                    linewidth=2,
                    ax=ax)
        ax.set_ylim(0, 1.05)
        ax.tick_params(axis='x', labelrotation=30)
        ax.set_xlabel("")
        ax.set_ylabel("Log$_{2}$(Abundance enrichment)")
        plt.tight_layout()
        plt.savefig(output[0])
