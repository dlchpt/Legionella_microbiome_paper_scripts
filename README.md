## Bacterial communities of premise plumbing systems in four European cities, and their association with culturable Legionella

#### Maria Scaturro, Federica Del Chierico, Yair Motro, Angeliki Chaldoupi, Anastasia Flountzi, Jacob Moran-Gilad, Antonietta Girolamo, Thomai Koutsiomani, Bozena Krogulska, Diane S. Lindsay, Renata Matuszewska, Georgios Papageorgiou, Katarzyna Pancer, Nikolaos Panoussis, Maria Cristina Rota, Soren Uldum, Emmanuel Velonakis, Dominique L. Chaput, Maria Luisa Ricci, ESCMID Study Group for Legionella Infections (ESGLI)

### Abstract
*Legionella* species are Gram negative, facultative, intracellular bacteria ubiquitous in natural and engineered water systems. Understanding the bacterial interactions underlying *Legionella*’s success in aquatic environments could be beneficial for control. We aimed to profile, by 16S rRNA amplicon sequencing, the bacterial communities in premise plumbing systems of buildings in four European cities (Copenhagen, Warsaw, Rome, Athens), and identify positive and negative associations of specific community members to culturable *Legionella*. Coarse taxonomic composition was similar across the four cities, but Copenhagen and Warsaw had richer, more diverse communities than Athens and Rome, with a greater number of city-specific ASVs. The cities had significantly different bacterial communities at the ASV level, with relatively few shared ASVs. Out of 5128 ASVs, 73 were classified as *Legionella*, and one or more of these were detected in most samples from each city (88.1% overall). Relative abundance of *Legionella* ASVs did not correlate with *Legionella* culture status. Overall, 44.2% of samples were *Legionella* culture positive: 71.4% in Warsaw, 62.2% in Athens, 22.2% in Rome, and 15.2% in Copenhagen. 54 specific ASVs and 42 genera had significant positive or negative associations with culturable *Legionella*. Negative associations included *Staphylococcus*, *Pseudomonas*, and *Acinetobacter*. Positive associations included several *Nitrospira* ASVs and one classified as *Nitrosomodaceae* oc32, ASVs in the amoeba-associated genera *Craurococcus-Caldovatus* and *Reyranella*, and the predatory genus *Bdellovibrio*. Some of these associations are well supported by laboratory studies, but others are the opposite of what was expected. This highlights the difficulties in translating pure culture results to complex real-world scenarios. However, these positive and negative associations held across the four cities, across multiple buildings and plumbing compartments. This is important because developing better control measures, including probiotic approaches, will require an understanding of ecological relationships that can be generalised across different engineered water systems.

### Importance
This study provides a snapshot of the diversity of microbial communities among premise plumbing systems in four European cities, providing new information on bacterial ASVs and genera that have positive or negative associations with culturable *Legionella* across a broad geographical range. This could inform studies aimed at confirming both in vitro and in real-world scenarios the role of other microbial community members in modulating *Legionella* proliferation, and will help in the development of probiotic approaches to controlling this opportunistic pathogen.


This repository details the analysis of the 16S rRNA amplicon sequence data.

Folder organisation follows the Schloss lab template, available from https://github.com/SchlossLab/new_project/releases/tag/0.14.


### Overview

	project
	|- README          # the top level description of content (this doc)
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- manuscript.Rmd    # executable Rmarkdown for this study, if applicable
	| |- references.bib # BibTeX formatted references
	| |- mbio.csl     # csl file to format references in the style of mBio
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- process/     # cleaned data, will not be altered once created;
	|
	|- code/          # any programmatic code to process data and make figures
	|
	|- results        # output from workflows and analyses
	| |- tables/      # text version of tables
	| |- figures/     # graphs, designated for manuscript figures
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary and additional analyses



### How to regenerate this repository
The metadata file is available here, `/data/raw/Metadata_main.txt`.

Raw sequencing reads are available from the European Nucleotide Archive under BioProject PRJEB52062. In order to run the analysis scripts from the start, these .fastq.gz files must be downloaded and placed in the three sequencing batch folders in `/data/raw/16S/`. All file names and their sequencing batch are listed in the metadata file. *I intend to write a bash script at some point that automates the ENA download and parsing into batch folders...*

Processed data, including the final cleaned data sets needed to run the analyses in `manuscript.Rmd` and generate all the figures (without needing to start from the raw .fastq files), are available on Figshare (https://doi.org/10.6084/m9.figshare.20685301.v1 for the ASV table, https://doi.org/10.6084/m9.figshare.20685370.v1 for the breakaway and ALDEx2 outputs). These should be placed in `/data/process/`. If the earlier steps are run (DADA2, etc.), their outputs will overwrite these Figshare files. 

Reference files used in this workflow were the latest Silva NR database release available at the time of analysis, version 138. There are two DADA2-formatted files:

* Main database, for overall classification, available at https://zenodo.org/record/3731176/files/silva_nr_v138_train_set.fa.gz?download=1
* Species assignment file, to refine classification to species level where possible, available at https://zenodo.org/record/3731176/files/silva_species_assignment_v138.fa.gz?download=1

These should be placed in the folder `/data/references/`.

#### Dependencies and locations
* R (v. 4.1.1) should be located in the user's PATH
* Main R packages (not exhaustive):
  * `ALDEx2`
  * `Biostrings`
  * `breakaway`
  * `dada2`
  * `decontam`
  * `factoextra`
  * `flextable`
  * `glmmTMB`
  * `here`
  * `knitr`
  * `lme4`
  * `lmerTest`
  * `lsmeans`
  * `microbiome`
  * `phyloseq`
  * `rmarkdown`
  * `tidyverse`
  * `UpSetR`
  * `vegan`