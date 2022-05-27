SAMPLES = ["t4_1_hrp_2", "t4_2_h2o2_1", "t5_1_hrp_2", "t5_2_h2o2_2"]
FORMATS = ["csv", "png"]

rule all:
    input: 
        expand("../resources/{sample}.{ext}",
               sample=SAMPLES, ext=FORMATS)

rule run_analysis:
    input:
        dataset="../resources/t4t5_12_frac_mass-spec_dataset.csv",
        mappings="../resources/uniprot2flybase.tab",
        annotations="../resources/localization_dataset.csv"
    params:
        sample="{sample}"
    output:
        "{sample}.{ext}"
    shell:
        "scripts/differential_analysis_mass-spec.py -d {input.dataset} -m {input.mappings} -a {input.annotations} -l {params.sample}"
