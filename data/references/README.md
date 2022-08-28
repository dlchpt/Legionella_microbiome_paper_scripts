# Reference files
Reference files used in this workflow were the latest Silva NR database release available at the time of analysis, version 138. There are two DADA2-formatted files:

* Main database, for overall classification, available at https://zenodo.org/record/3731176/files/silva_nr_v138_train_set.fa.gz?download=1
* Species assignment file, to refine classification to species level where possible, available at https://zenodo.org/record/3731176/files/silva_species_assignment_v138.fa.gz?download=1

These should be placed in the folder `/data/references/`. Ensure they are named `silva_nr_v138_train_set.fa.gz` and `silva_species_assignment_v138.fa.gz`.