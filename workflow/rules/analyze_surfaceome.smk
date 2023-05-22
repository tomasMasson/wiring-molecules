from Bio import SeqIO
import pandas as pd


rule all:
    input: 
        "../results/surfaceome_analysis/dmel_extracellular_proteome.fasta",
        "../results/surfaceome_analysis/dmel_extracellular_domains.csv",
        "../results/surfaceome_analysis/dmel_adhesion_proteome.csv"


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


rule download_uniprot_data:
    params:
        proteome="http://ftp.flybase.net/releases/current/dmel_r6.51/fasta/dmel-all-translation-r6.51.fasta.gz",
        gene_ontology="http://ftp.flybase.net/releases/current/precomputed_files/go/gene_association.fb.gz",
        uniprot_mappings="http://ftp.flybase.net/releases/current/precomputed_files/genes/fbgn_NAseq_Uniprot_fb_2023_02.tsv.gz"
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
        "grep -e Immunoglobulin -e Leucine -e fn3 {input} | awk '{{print $4}}' | sort | uniq > {output}"


rule move_files_to_results:
    input:
        prot="dmel_extracellular_proteome.fasta",
        dom="dmel_extracellular_domains.csv",
        adh="dmel_adhesion_proteome.csv"
    output:
        prot="../results/surfaceome_analysis/dmel_extracellular_proteome.fasta",
        dom="../results/surfaceome_analysis/dmel_extracellular_domains.csv",
        adh="../results/surfaceome_analysis/dmel_adhesion_proteome.csv"
    shell:
        """
        mv {input.prot} {output.prot} && \
        mv {input.dom} {output.dom} && \
        mv {input.adh} {output.adh}
        """
