# Command line snippet to retrieve the protein sequences for a list of identifiers

`for i in (cat id.csv); rsync -av rsync://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/metazoa/fasta/$i/pep/ datasets/$i; end`
