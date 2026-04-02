# Supporting data analysis pipeline for the paper: "A Molecular Grammar to Assemble Motion Vision Circuits in *Drosophila*"

This repository contains the code and data required to reproduce the analysis performed in [INSERT PREPRINT LINK].

## Cell-surface proteomics profiling in T4 and T5 neurons

This analysis in consolidated on the snakefile `~/wiring-molecules/workflow/rules/ratiometric_analysis.smk`. After performing mass-spec runs on our cell-surface enriched samples, raw data from the instrument was deconvolved as described in the paper methods and stored as a table (`/wiring-molecules/resources/t4vt5_pl_dataset.csv`) where protein hits are on the rows and each channel from the TMT quantification is on a different column. From this point, we performed a ratiometric analysis, as described in "Cell-Surface Proteomic Profiling in the Fly Brain Uncovers Wiring Regulators" PMID:31955847, using the script `/wiring-molecules/workflow/scripts/ratiometric_classifier.py`. Briefly, GO terms from FlyBase (`/wiring-molecules/resources/gene_association.fb`) and extracellular proteins from FlyXCDB (`/wiring-molecules/resources/flyxcdb_data.csv`) are collated as gold standard and used to score the False Discovery Rate (False Positives / (True Positives + False Positives)). For the final subset, we keep all proteins with a signal ratio (intensity / control intensity) bigger than the ratio of the data point where the difference between true positive rate and false positive rate is maximum; the intersection between sample replicates give us the final set (`/wiring-molecules/results/t4t5_surface_proteomics/t4t5_consensus_surfaceome.csv`). For each comparison against the four control conditions, we plot the Receiver Operating Characteristic (ROC) curve, the distribution of intensity ratios for the true positives and false positives and the difference between true positive rate and false positive rate as function of the intensity ratio (results are stored in `/wiring-molecules/results/t4t5_surface_proteomics`).

After ratiometric analysis, correlation plots between replicates (`/wiring-molecules/results/t4t5_surface_proteomics/replicates_pearson_correlation.svg`) and a PCA (`/wiring-molecules/results/t4t5_surface_proteomics/pca_analysis.svg`) are performed to check the quality of the raw data. Scripts for these are stored inside `wiring-molecules/workflow/scripts`.



## Summary

Instinctive behaviors depend on neuronal circuits that are hardwired into animal genomic information. How this molecular process is achieved remains obscure, but a necessary ingredient to achieve these precise wiring is the presence of recognition codes at the neuronal surface. Although multiple families of cell adhesion molecules have been implicated in synaptic specificity, a systematic picture of how these repertoires instruct partner matching remains vacant. Here, we aim to uncover the molecular assembly principles that support instinctive behaviors. Delineating these rules will open the door to program novel neuronal circuits and behaviors through synthetic approaches.

## Uncovering the Molecular Grammar of Brain Wiring in *Drosophila*

### *In silico* prediction of the neuronal adhesion surfaceome

The workflow `workflow/rules/analyze_drosophila_surfaceome.smk` uses data from DeepLoc2 localization prediction, Gene Ontology (GO) terms classification and the FlyXCDB to consolidate a comprehensive set of surface proteins in _Drosophila_, the so called Surfaceome. DeepLoc2, GO and FlyXCDB datasets are provided as precomputed files in the resources folder (`resources/dmel_proteome_localization_deeploc2.csv`, `resources/gene_association.fb`, `resources/flyxcdb_data.csv`). UniProt → FlyBase mapping is based on `resources/fbgn_NAseq_Uniprot_fb_2024_06.tsv`.

After selecting protein belonging to the plasma membrane, scRNAseq data (`resources/FlyCellAtlas_slimmed_gene_expression_fb_2024_06.tsv`) is used to keep only genes that are expressed in either neuron or interneuron classes (cluster fraction expressing the gene > 0). Then, `db="resources/Pfam-A.hmm` database is queried to perform domain annotation through hmmscan.

`snakemake -s workflow/rules/analyze_drosophila_surfaceome.smk --cores 6`

Outputs include a multifasta with all the surface proteins in _Drosophila (`results/synapse_interaction_maps/drosophila_surfaceome.faa`), a surfaceome subset that is expressed in neurons/interneurons (`results/synapse_interaction_maps/drosophila_neuronal_surfaceome.faa`) domain hits for those proteins (`results/synapse_interaction_maps/drosophila_neuronal_surfaceome_domains_search.csv`) and a list of proteins selected as the wiring surfaceome based on the presence Immunoglobulin/Leucine/Fibronectin 3/EGF domains (`results/synapse_interaction_maps/drosophila_wiring_surfaceome.csv`).

# Adaptation of adhesion molecules involved in neuronal circuits wiring

## Estimating cotransitions across wiring molecules in Diptera

### Usage

To run this part of the pipeline you'll need to move into the `owrkflow` and run the Snakemake rule `rules/get_orthodb_genecounts_matrix.smk` together with the `rules/config.yaml` configuration file:

`snakemake --snakefile rules/get_orthodb_cotransitions.smk --configfile rules/config.yaml --cores 4`

This pipeline will download all the data required from `https://data.orthodb.org/download/` and FlyBase, using them to compute the contransition scores for each orthogroup annotated at the Diptera level (NCBI taxon ID 7147). After completing the run, results should be presents inside `results/orthodb_cotransitions/`.
