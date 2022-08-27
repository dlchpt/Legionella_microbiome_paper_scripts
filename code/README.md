# Analysis scripts for amplicon data
All the R scripts used to go from raw data (.fastq.gz files) to final figures and results will be stored here. Discrete steps in data analysis will be kept in separate scripts to make it easier to re-run and/or change settings, and I will eventually make a batch file that calls on each of these scripts sequentially so that the whole analysis can be run in a single command.

The main steps are:
1. `DADA2_16S.R`, which goes from raw data to an amplicon sequence variant (ASV) table.
1. `prune_filter_final_dataset.R`, which removes failed libraries and non-target sequences, applies a prevalence filter to remove likely sequencing errors, and assigns species where possible.
1. `decontam_16S.R`, which removes likely contaminants introduced during DNA extraction and library prep.
1. `generate_tree_16S.R` takes the decontaminated data and builds a phylogenetic tree of the ASV sequences. This step is computationally intensive and must be run on a cluster.
1. With the clean streamlined ASV tables, I can then move onto proper analysis: alpha/beta diversity, differential abundance of taxa, core microbiome analysis, modelling with phys/chem parameters, etc.