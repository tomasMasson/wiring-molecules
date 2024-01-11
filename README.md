# Adaptation of adhesion molecules involved in neuronal circuits wiring

## Estimating cotransitions across wiring molecules in Diptera

### Usage

To run this part of the pipeline you'll need to move into the `owrkflow` and run the Snakemake rule `rules/get_orthodb_genecounts_matrix.smk` together with the `rules/config.yaml` configuration file:

`snakemake --snakefile rules/get_orthodb_cotransitions.smk --configfile rules/config.yaml --cores 4`

This pipeline will download all the data required from `https://data.orthodb.org/download/` and FlyBase, using them to compute the contransition scores for each orthogroup annotated at the Diptera level (NCBI taxon ID 7147). After completing the run, results should be presents inside `results/orthodb_cotransitions/`.
