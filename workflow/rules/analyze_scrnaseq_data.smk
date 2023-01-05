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
        expand("{sample}.png", sample=SAMPLES),
        expand("{sample}_no_expr_t4.csv", sample=SAMPLES),
        expand("{sample}_no_expr_t5.csv", sample=SAMPLES),

rule run_analysis:
    input:
        folder="t4t5_scrnaseq_data/",
        t4_proteome="workflow/t4_consensus_surfaceome.csv",
        t5_proteome="workflow/t5_consensus_surfaceome.csv",
    params:
        sample="{sample}",
    output:
        "{sample}.png",
        "{sample}_no_expr_t4.csv",
        "{sample}_no_expr_t5.csv"
    shell:
        """
        ./process_scrnaseq.py --folder {input.folder} --sample {params} --t4_proteome {input.t4_proteome} --t5_proteome {input.t5_proteome}
        """
