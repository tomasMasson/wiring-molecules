import scanpy as sc
import numpy as np
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


### Auxiliary functions

def process_scrnaseq_dataset(folder: str, prefix: str, annotation: str):
    # adata = sc.read_10x_mtx("zipursky_t4t5_scrnaseq_dataset", var_names="gene_ids", prefix="GSM3592259_T4T5_24_")
    # Ingest data
    adata = sc.read_10x_mtx(folder,
                            var_names="gene_ids",
                            prefix=prefix)
    # Discard cells with low gene counts
    sc.pp.filter_cells(adata, min_genes=200)
    # Discard genes with low cell counts
    sc.pp.filter_genes(adata, min_cells=3)
    # Normalize counts (10000 factor)
    sc.pp.normalize_total(adata, target_sum=1e4)
    # annotation = pd.read_csv("zipursky_t4t5_scrnaseq_dataset/GSM3592259_clust24.txt.gz", sep="\t")
    # Annotate neuron types
    annotation = pd.read_csv(annotation, sep="\t")
    adata.obs["Celltype"] = adata.obs.index.str[:-2]
    adata.obs.replace(dict(zip(annotation["barcode"], annotation["cluster"])), inplace=True)
    # Transform matrix into a DataFrame
    df = pd.DataFrame()
    for row, idx in adata.obs.groupby(["Celltype"]):
        if row == "filtered":
            pass
        else:
            fadata = adata[idx.index, :].to_df()
            means = fadata.mean(axis=0)
            celltypes = [row] * len(means)
            dff = pd.DataFrame({"FBgn": means.index, "Celltype": celltypes, "Expression": means})
            df = pd.concat([df, dff])
    return df

### End auxiliary functions


rule all:
    input:
        "../results/t4t5_surface_proteomics/t4t5_consensus_surfaceome.csv",
        "../results/t4t5_surface_proteomics/pca_analysis.svg",
        "../results/t4t5_surface_proteomics/venn_diagram.svg",
        "../results/t4t5_surface_proteomics/t4t5_surfaceome_physical_interactions.csv",
        "../results/t4t5_surface_proteomics/cell_adhesion_molecules_abundance.svg",
        "../results/t4t5_surface_proteomics/cell_adhesion_molecules_expression.svg"


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
        fig, ax = plt.subplots(figsize=(4,4))
        dff = df.replace({True: "T4 Wiring",
                          False: "T4 Dataset"})
        sns.boxplot(x="CAM",
                    y="T4_signal",
                    data=dff[dff["T4_signal"] != 0],
                    color="#f1a340ff",
                    showcaps=False,
                    linewidth=2,
                    ax=ax)
        dff = df.replace({True: "T5 Wiring",
                          False: "T5 Dataset"})
        sns.boxplot(x="CAM",
                    y="T5_signal",
                    data=dff[dff["T5_signal"] != 0],
                    color="#998ec3ff",
                    showcaps=False,
                    linewidth=2,
                    ax=ax)
        ax.set_ylim(0, 1.05)
        ax.tick_params(axis='x', labelrotation=15)
        ax.set_xlabel("")
        ax.set_ylabel("Log$_{2}$(Abundance enrichment)")
        plt.tight_layout()
        plt.savefig(output[0])


rule plot_non_expressed_adhesion_molecules:
    input:
        "../results/t4t5_surface_proteomics/t4t5_consensus_surfaceome.csv",
        "../resources/expanded_cams_list.csv"
    params:
        dataset_24h="../resources/zipursky_t4t5_24h_scrnaseq_dataset/",
        prefix_24h="GSM3592259_T4T5_24_",
        clusters_24h="../resources/zipursky_t4t5_24h_scrnaseq_dataset/GSM3592259_clust24.txt.gz",
        dataset_48h="../resources/zipursky_t4t5_48h_scrnaseq_dataset/",
        prefix_48h="GSM3592260_T4T5_48_",
        clusters_48h="../resources/zipursky_t4t5_48h_scrnaseq_dataset/GSM3592260_clust48.txt.gz"
    output:
        "../results/t4t5_surface_proteomics/cell_adhesion_molecules_expression.svg",
        "../results/t4t5_surface_proteomics/t4t5_expressed_cell_adhesion_molecules_24h.csv",
        "../results/t4t5_surface_proteomics/t4t5_expressed_cell_adhesion_molecules_48h.csv"
    run:
        # Load proteomic data
        df = pd.read_csv(input[0])
        # Load cell adhesion molecules data
        cam = pd.read_csv(input[1])
        # Process 24h expression data
        exp_24h = process_scrnaseq_dataset(params["dataset_24h"], params["prefix_24h"], params["clusters_24h"])
        cams_24 = set(exp_24h[(exp_24h["FBgn"].isin(cam["FBgn"]))&(exp_24h["Expression"] > 1)]["FBgn"])
        with open(output[1], "w") as fh:
            for i in cams_24:
                fh.write(f"{i}\n")
        exp_24h = exp_24h[exp_24h["FBgn"].isin(df["Protein"])]
        exp_24h = exp_24h[exp_24h["FBgn"].isin(cam["FBgn"])]
        # Process 48h expression data
        exp_48h = process_scrnaseq_dataset(params["dataset_48h"], params["prefix_48h"], params["clusters_48h"])
        cams_48 = set(exp_48h[(exp_48h["FBgn"].isin(cam["FBgn"]))&(exp_48h["Expression"] > 1)]["FBgn"])
        with open(output[2], "w") as fh:
            for i in cams_48:
                fh.write(f"{i}\n")
        cams_48.intersection(set(df["Protein"]))
        exp_48h = exp_48h[exp_48h["FBgn"].isin(df["Protein"])]
        exp_48h = exp_48h[exp_48h["FBgn"].isin(cam["FBgn"])]
        # Merge 24h and 48h datasets
        df = pd.concat([exp_24h, exp_48h])
        df = df.groupby("FBgn")["Expression"].max()
        df = df.sort_values(ascending=False)
        # Plot expression data
        sns.set(font_scale=1.2)
        sns.set_style("white")
        fig, ax = plt.subplots(figsize=(4,4))
        sns.scatterplot(x=range(len(df)), y=np.log1p(df.values), color="gray", linewidth=0, s=16, ax=ax)
        ax.set_ylim(-.2, 4.6)
        ax.set_xlabel("Ranked Cell Adhesion Molecules")
        ax.set_ylabel("Log1p(Maximum Expression)")
        plt.tight_layout()
        plt.axhline(.5, color="black", linestyle="--")
        plt.savefig(output[0])
