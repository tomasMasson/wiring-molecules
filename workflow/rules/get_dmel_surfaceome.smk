import pandas as pd
from Bio import SeqIO
from Bio import SearchIO
from itertools import combinations


rule all:
    input: 
        "../results/surfaceome_analysis/dmel_extracellular_proteome.fasta",
        "../results/surfaceome_analysis/dmel_extracellular_domains.csv",
        "../results/surfaceome_analysis/dmel_extracellular_architecture.csv",
        "../results/surfaceome_analysis/dmel_adhesion_proteome.csv",
        "../results/surfaceome_analysis/dmel_t4t5_candidate_adhesion_interactions.fasta"

# Start auxiliar functions

def extract_protein_subset(protein_list: list, proteome: str, outfile: str) -> str:
    "Returns a protein subset from a multifasta"
    # Check if the output provided
    if outfile is None:
        outfile = f"dmel_membrane_proteins.fasta"
    # Store protein headers after processing
    idlist = []
    with open(outfile, "w") as fh:
        for seq in SeqIO.parse(proteome, "fasta"):
            # Get FlyBase (FBgn) identifier
            header = seq.description.split("parent=")[1].split(",")[0]
            # Check if already processed
            if header in protein_list:
                if header not in idlist:
                    # Write output
                    fh.write(f">{header}\n{seq.seq}\n")
                    idlist.append(header)
    return outfile

# End auxiliar functions


rule download_uniprot_data:
    params:
        proteome="http://ftp.flybase.net/releases/FB2023_04/dmel_r6.53/fasta/dmel-all-translation-r6.53.fasta.gz",
        gene_ontology="http://ftp.flybase.net/releases/FB2023_04/precomputed_files/go/gene_association.fb.gz",
        uniprot_mappings="http://ftp.flybase.net/releases/FB2023_04/precomputed_files/genes/fbgn_NAseq_Uniprot_fb_2023_04.tsv.gz"
    output:
        temp("dmel_proteome.fasta"),
        temp("gene_association.fb"),
        temp("uniprot_mappings.tsv")
    shell:
        """
        wget {params.proteome} -O dmel_proteome.fasta.gz && \
        gunzip dmel_proteome.fasta.gz && \
        wget {params.gene_ontology} -O gene_association.fb.gz &&\
        gunzip gene_association.fb.gz
        wget {params.uniprot_mappings} -O uniprot_mappings.tsv.gz &&\
        gunzip uniprot_mappings.tsv.gz
        """


rule filter_deeploc2_extracellular_proteins:
    input:
        "../resources/dmel_proteome_localization_deeploc2.csv",
    output:
        temp("deeploc2_extracellular.csv")
    shell:
        "grep -e 'Cell membrane' -e 'Extracellular' {input} > {output}"


rule rename_deeploc2_predictions:
    input:
        deepl="deeploc2_extracellular.csv",
        fbgn="uniprot_mappings.tsv"
    output:
        temp("deeploc2_extracellular_fbgn.csv")
    run:
        df = pd.read_csv(input.deepl,
                         usecols=[0, 1],
                         )
        df.Protein_ID = df.Protein_ID.str.split("|", expand=True)[1]
        labels = pd.read_csv(input.fbgn, sep="\t", names=["FBgn", "UniProt"], skiprows=5, skipfooter=1, usecols=[2, 5]).dropna()
        labels = labels[labels.UniProt.isin(df.Protein_ID)]
        df.replace(dict(zip(labels.UniProt, labels.FBgn)), inplace=True)
        df.to_csv(output[0], index=False)


rule extract_extracellular_proteins:
    input:
        proteome="dmel_proteome.fasta",
        gene_ontology="gene_association.fb",
        deeploc="deeploc2_extracellular_fbgn.csv",
    output:
        "dmel_extracellular_proteome.fasta"
    run:
        df = pd.read_csv(input.gene_ontology,
                         sep="\t",
                         names=["FBgn", "Gene", "GO"],
                         usecols=[1, 2, 4],
                         skiprows=5)
        ex_proteome = []
        # Plasma membrane GO
        ex_proteome += list(df[df.GO == "GO:0005886"].FBgn)
        # Cell-cell junction GO
        ex_proteome += list(df[df.GO == "GO:0005911"].FBgn)
        # Extracellular region GO
        ex_proteome += list(df[df.GO == "GO:0005576"].FBgn)
        # Deeploc2 prediction
        deeploc = pd.read_csv(input.deeploc, names=["Protein_ID", "Localizations"])
        ex_proteome += list(deeploc.Protein_ID)
        extract_protein_subset(ex_proteome, input.proteome, output[0])


rule search_pfam_domains:
    input:
        "dmel_extracellular_proteome.fasta"
    output:
        "dmel_extracellular_domains.csv"
    params:
        db="../resources/Pfam-A.hmm"
    shell:
        """
        hmmscan --domtblout {output} --noali --cpu 4 -E 0.001 --domE 0.001 {params.db} {input}
        """


rule extract_adhesion_proteome:
    input:
        "dmel_extracellular_domains.csv"
    output:
        "dmel_adhesion_proteome.csv"
    shell:
        "grep -e Immunoglobulin -e Leucine -e fn3 -e EGF {input} | awk '{{print $4}}' | sort | uniq > {output}"


rule extract_protein_architecture:
    input:
        "dmel_extracellular_domains.csv"
    output:
        "dmel_extracellular_architecture.csv"
    run:
        domains = {}
        search = SearchIO.parse(input[0],
                                "hmmsearch3-domtab")
        for query in search:
            dom = []
            for hit in query:
                dom.append(hit.id)
            dom = ("/".join(dom))
            domains[query.id] = dom
        df = pd.DataFrame.from_dict(domains, orient="index", columns=["Domains"])
        df.to_csv(output[0])


rule extract_t4t5_candidate_adhesion_interactions:
    input:
        adhesome="dmel_adhesion_proteome.csv",
        t4t5="../results/pl_proteomics/pl_results/t4t5_consensus_surfaceome.csv",
        fastas="dmel_extracellular_proteome.fasta"
    output: "dmel_t4t5_candidate_adhesion_interactions.fasta"
    run:
        adhesome = pd.read_csv(input.adhesome, names=["Gene"])
        t4t5 = pd.read_csv(input.t4t5)
        interactions = list(combinations(t4t5[t4t5.Protein.isin(adhesome.Gene)].Protein, 2))
        seqs = SeqIO.to_dict(SeqIO.parse(input.fastas, "fasta"))
        with open(output[0], "w") as fh:
            for interaction in interactions:
                seq1 = seqs[interaction[0]]
                seq2 = seqs[interaction[1]]
                fh.write(f">{seq1.id}_{seq2.id}\n{seq1.seq}:{seq2.seq}\n")


rule move_files_to_results:
    input:
        prot="dmel_extracellular_proteome.fasta",
        dom="dmel_extracellular_domains.csv",
        arc="dmel_extracellular_architecture.csv",
        adh="dmel_adhesion_proteome.csv",
        ints="dmel_t4t5_candidate_adhesion_interactions.fasta"
    output:
        prot="../results/surfaceome_analysis/dmel_extracellular_proteome.fasta",
        dom="../results/surfaceome_analysis/dmel_extracellular_domains.csv",
        arc="../results/surfaceome_analysis/dmel_extracellular_architecture.csv",
        adh="../results/surfaceome_analysis/dmel_adhesion_proteome.csv",
        ints="../results/surfaceome_analysis/dmel_t4t5_candidate_adhesion_interactions.fasta"
    shell:
        """
        mv {input.prot} {output.prot} && \
        mv {input.dom} {output.dom} && \
        mv {input.arc} {output.arc} && \
        mv {input.adh} {output.adh} && \
        mv {input.ints} {output.ints}
        """
