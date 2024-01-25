import numpy as np
import pandas as pd
from ete3 import NCBITaxa, Tree

# Start auxiliar functions



# End auxiliar functions

TAXONS = config["TAXON_ID"]


rule all:
    input:
        "../results/orthodb_cotransitions/orthodb_og2genes.tab.gz",
        "../results/orthodb_cotransitions/orthodb_gene_annotations.tab.gz",
        "../results/orthodb_cotransitions/flybase_annotations.tsv.gz",
        expand("../results/orthodb_cotransitions/{taxa}_specific_genecounts.csv", taxa=TAXONS),
        expand("../results/orthodb_cotransitions/{taxa}_specific.treefile", taxa=TAXONS),
        expand("../results/orthodb_cotransitions/{taxa}_specific_cotransitions.tsv", taxa=TAXONS)


rule download_orthodb_orthogroups:
    params: config["ORTHODB_OG2GENES"]
    output: "orthodb_og2genes.tab.gz"
    shell:  "wget {params} -O {output}"


rule download_orthodb_gene_annotations:
    params: config["ORTHODB_GENES_ANN"]
    output: "orthodb_gene_annotations.tab.gz"
    shell:  "wget {params} -O {output}"


rule download_FlyBase_annotations:
    params: config["ANNOTATION_DB"]
    output: "flybase_annotations.tsv.gz"
    shell:  "wget {params} -O {output}"


rule parse_taxon_orthogroups:
    input:  "orthodb_og2genes.tab.gz"
    params: taxa="{taxa}"
    output: temp("{taxa}_specific_orthogroups.tab")
    shell:  "zcat {input} | grep 'at{params.taxa}' > {output}"


rule parse_reference_orthodb_mappings:
    input:  "orthodb_og2genes.tab.gz"
    params:
        reference=config["REFERENCE_SPECIE"],
        taxa="{taxa}"
    output: temp("{taxa}_reference_og2genes_mapping.tab")
    shell:
        """
        zcat {input} | \
        grep -P '\t{params.reference}_0:' | \
        grep 'at{params.taxa}' > {output}
        """


rule parse_orthodb_gene_annotations:
    input:  "orthodb_gene_annotations.tab.gz"
    params: reference=config["REFERENCE_SPECIE"]
    output: temp("orthodb_reference_gene_annotations.tab")
    shell:
        """
        zcat {input} | \
        grep -w '{params.reference}_0' > {output}
        """


rule get_orthodb_genecounts_matrix:
    input:
        og="{taxa}_specific_orthogroups.tab",
        reference="{taxa}_reference_og2genes_mapping.tab",
        uniprot="orthodb_reference_gene_annotations.tab",
        flybase="flybase_annotations.tsv.gz"
    output: "{taxa}_specific_genecounts.csv"
    run:
        # Load orthogroups data
        data = pd.read_csv(input.og,
                           names=["og", "gene"],
                           sep="\t")
        # Extract NCBI's species names
        data.loc[:, "species"] = data.gene.str.split("_0", expand=True)[0]
        # Compute gene counts across orthogroups/species
        genecounts = data.groupby(["og", "species"]).size().unstack(fill_value=0)

        # Load reference orthogroup mapping
        ref_data = pd.read_csv(input.reference,
                               names=["og", "gene"],
                               sep="\t")
        # Load reference orthogroup annotations
        ref_ann = pd.read_csv(input.uniprot,
                              usecols=[0, 4],
                              names=["gene", "uniprot"],
                              sep="\t")
        # Load reference Flybase annotations
        ref_flybase = pd.read_csv(input.flybase,
                                  skiprows=5,
                                  usecols=[2, 5],
                                  names=["fbgn", "uniprot"],
                                  sep="\t")

        # Rename OrthoDB OGs -> Gene
        genecounts.rename(dict(zip(ref_data.og, ref_data.gene)), axis="index", inplace=True)
        # Rename OrthoDB Gene -> UniProt ID
        genecounts.rename(dict(zip(ref_ann.gene, ref_ann.uniprot)), axis="index", inplace=True)
        # Rename UniProt ID -> FlyBase ID
        genecounts.rename(dict(zip(ref_flybase.uniprot, ref_flybase.fbgn)), axis="index", inplace=True)
        if "180454" in genecounts.columns:
            # Remove problematic leaf (wrongly parsed by ete)
            genecounts.drop("180454", axis=1, inplace=True)
        if "278856" in genecounts.columns:
            # Remove problematic leaf (wrongly parsed by ete)
            genecounts.drop("278856", axis=1, inplace=True)
        if "121224" in genecounts.columns:
            # Remove problematic leaf (wrongly parsed by ete)
            genecounts.drop("121224", axis=1, inplace=True)
        if "441943" in genecounts.columns:
            # Rename problematic leaf (441943 -> 2961670 N. virginianus)
            genecounts.rename({"441943": 2961670}, axis="columns", inplace=True)
        # Save dataframe
        genecounts.to_csv(output[0])


rule extract_ncbi_taxonomic_constraints:
    input:  "{taxa}_specific_genecounts.csv"
    output: temp("{taxa}_specific_ncbi_taxonomy.nwk")
    run:
        genecounts = pd.read_csv(input[0], index_col="og")
        ncbi = NCBITaxa()
        taxa = genecounts.columns
        taxa = [t.split("_")[0]
                for t in taxa]
        tree = ncbi.get_topology(taxa)
        leaves = [leaf.get_leaf_names()[0]
                  for leaf in tree.get_leaves()]
        print(len(taxa))
        print(len(leaves))
        tree.write(outfile=output[0], format=9)


rule generate_genecounts_fasta:
    input:  "{taxa}_specific_genecounts.csv"
    output: temp("{taxa}_specific_profiles.fasta")
    run:
        data = pd.read_csv(input[0], index_col="og")
        data = data[data.sum(axis=1) > 4]
        with open(output[0], "w") as fh:
            for row in data.T.iterrows():
                species = str(row[0]).split("_")[0]
                profile = pd.get_dummies(row[1] > 0).astype(int).iloc[:,1].astype("str")
                profile = "".join(list(profile))
                fh.write(f">{species}\n{profile}\n")


rule infer_ml_phylogeny:
    input:  "{taxa}_specific_profiles.fasta",
            "{taxa}_specific_ncbi_taxonomy.nwk"
    params: "{taxa}_specific"
    output: "{taxa}_specific.treefile"
    shell:
        """
        iqtree -redo -s {input[0]} -g {input[1]} -m GTR2+FO+R5 --prefix {params[0]}
        """


rule compute_gene_cotransition:
    input:  "{taxa}_specific_genecounts.csv",
            "{taxa}_specific.treefile"
    output: "{taxa}_specific_cotransitions.tsv"
    run:
        # Load treefile and extract leaves order
        tree = Tree(input[1])
        tree.ladderize(direction=0)
        leaves_order = tree.get_leaf_names()
        # Load genecount matrix
        data = pd.read_csv(input[0], index_col="og")
        data = data[data.sum(axis=1) > 4]
        data = data[leaves_order]
        ngenes, norgs = len(data.index), len(data.columns)
        # Compute cotransitions (as Dembech et. al)
        transitions = data.diff(axis=1).applymap(lambda x: 1 if x>1 else -1 if x <-1 else x).values.tolist()
        # sets for fast comparison
        t01 = [set(np.nonzero(row > 0)[0]) for row in np.array(transitions)] # 0->1 transitions
        t10 = [set(np.nonzero(row < 0)[0]) for row in np.array(transitions)] # 1->0 transitions
        tt = [len(a | b) for a,b in zip(t01,t10)] #total transitions

        with open(output[0], "w") as fh:
            fh.write(f"Orthogroup1\tOrthogroup2\torgs\tt1\tt2\tc\td\tk\n")
            counter = 0
            for i in range(ngenes-1):
                for j in range(i+1,ngenes):
                    concordant = len(t01[i] & t01[j]) + len(t10[i] & t10[j])
                    discordant = len(t10[i] & t01[j]) + len(t01[i] & t10[j])
                    k = concordant - discordant
                    if abs(k) >= 0:
                        fh.write(f"{data.index[i]}\t{data.index[j]}\t{norgs}\t{tt[i]}\t{tt[j]}\t{concordant}\t{discordant}\t{k}\n")
                        counter += 1


rule copy_results:
    input:
        "orthodb_og2genes.tab.gz",
        "orthodb_gene_annotations.tab.gz",
        "flybase_annotations.tsv.gz",
        expand("{taxa}_specific_genecounts.csv", taxa=TAXONS),
        expand("{taxa}_specific.treefile", taxa=TAXONS),
        expand("{taxa}_specific_cotransitions.tsv", taxa=TAXONS)
    output:
        "../results/orthodb_cotransitions/orthodb_og2genes.tab.gz",
        "../results/orthodb_cotransitions/orthodb_gene_annotations.tab.gz",
        "../results/orthodb_cotransitions/flybase_annotations.tsv.gz",
        expand("../results/orthodb_cotransitions/{taxa}_specific_genecounts.csv", taxa=TAXONS),
        expand("../results/orthodb_cotransitions/{taxa}_specific.treefile", taxa=TAXONS),
        expand("../results/orthodb_cotransitions/{taxa}_specific_cotransitions.tsv", taxa=TAXONS)
    shell: "mv {input} ../results/orthodb_cotransitions/ && rm *.gz *.iqtree *.log *.parstree"
