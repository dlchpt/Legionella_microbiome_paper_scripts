# Further analyses and manuscript preparation
This folder contains the main Rmarkdown file (`manuscript.Rmd`) used for further analyses and to generate the draft manuscript.

It starts with the phyloseq object `phyloseq_after_pruning_filtering.rds` generated at the end of the script `02_prune_filter_final_dataset.R`. If previous scripts haven't been run from the raw sequence files, the necessary phyloseq object can be downloaded from Figshare (https://doi.org/10.6084/m9.figshare.20685301.v1) and should be placed in the folder `/data/process/`. 

Some steps in this Rmd file are time-consuming, notably the estimation of population richness with `breakaway` and the differential abundance analysis with `ALDEx2`. To speed up the compilation of the manuscript, the data objects generated by these steps can instead be downloaded from Figshare (https://doi.org/10.6084/m9.figshare.20685370.v1) and should be placed in the folder `/data/process/`. If these downloaded files are being used instead of being computed (or if the Rmd file is being run again after these files have already been generated), the breakaway and ALDEx2 lines in the `manuscript.Rmd` file should be commented out (lines 215, 218, 814, 816, 861, 863). Otherwise, these will be computed again and will overwrite any downloaded or previously-computed files.