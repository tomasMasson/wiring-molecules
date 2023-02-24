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
        expand("{sample}_protein_expression.csv", sample=SAMPLES)

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
