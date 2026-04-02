DESPLAN_COUNTS = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142787/suppl/GSE142787%5FMixture%5Fmodeling%2Exlsx"
DESPLAN_LABELS = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2879-3/MediaObjects/41586_2020_2879_MOESM4_ESM.xlsx"
SAMPLES = ["GSM4451487_24h_APF_1",
           "GSM4451488_24h_APF_2",
           "GSM4451489_36h_APF_1",
           "GSM4451490_36h_APF_2",
           "GSM4451491_48h_APF_1",
           "GSM4451492_48h_APF_2",
           "GSM4451493_60h_APF_1",
           "GSM4451494_60h_APF_2",
           "GSM4451495_72h_APF_1",
           "GSM4451496_72h_APF_2",
           "adult"]


rule all:
    input: 
        "desplan_optic_lobe_atlas_counts.xlsx",
        "desplan_optic_lobe_atlas_labels.xlsx"


rule fetch_datasets:
    output: 
        counts="desplan_optic_lobe_atlas_counts.xlsx",
        labels="desplan_optic_lobe_atlas_labels.xlsx"
    shell:
        """
        wget {DESPLAN_COUNTS} -O {output.counts}
        wget {DESPLAN_LABELS} -O {output.labels}
        """

sheets = pd.read_excel("", sheet=None)
timepoints = sheets.keys()
t4/t5 ab = 234
t4/t5 cd = 235
labels = pd.read_csv("fbgn_annotation_ID.tsv", sep="\t", skiprows=4)
labs = dict(zip(labels["##gene_symbol"], labels["primary_FBgn#"]))



rule run_analysis:
    input:
        folder="../resources/t4t5_scrnaseq_data/",
        t4_proteome="../results/pl_proteomics/pl_results/t4_consensus_surfaceome.csv",
        t5_proteome="../results/pl_proteomics/pl_results/t5_consensus_surfaceome.csv"
    params:
        sample="{sample}",
    output:
        "{sample}_protein_expression.csv"
    shell:
        """
        ./scripts/process_scrnaseq.py --folder {input.folder} --sample {params} --t4_proteome {input.t4_proteome} --t5_proteome {input.t5_proteome}
        """
