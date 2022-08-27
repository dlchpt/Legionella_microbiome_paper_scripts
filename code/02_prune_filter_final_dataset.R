# Additional steps to prepare data for analysis
# * Remove failed library (sample_sums < 2000)
# * Prevalence filter to remove rarest ASVs (more likely to be sequencing errors)
# * Remove non-target sequences (mitochondria, chloroplasts, eukaryotes)
# * Try species taxonomy step once the previous two steps reduced the number of ASVs



# Load 16S V4 data
V4 <- readRDS(here("data", "process", "ps_object_after_dada2.rds"))

# Remove libraries smaller than 2000 reads
V4.trim <- prune_samples(sample_sums(V4) >= 2000, V4)
V4.trim <- prune_taxa(taxa_sums(V4.trim) > 0, V4.trim)

# Prevalence filter, keeping only ASVs present in 3 or more samples 
V4.filt <- filter_taxa(V4.trim, function(x) sum(x > 0) > 2, TRUE)

# Keep only taxa that are classified as Bacteria at the Domain level
V4.domain <- subset_taxa(V4.filt, Kingdom == "Bacteria")
V4.domain # 5140 taxa and 176 samples

# Remove mitochondria and chloroplasts
# subset_taxa defaults to removing NAs too, so need to specify that they be kept
V4.no.mito <- subset_taxa(V4.domain, (Family != "Mitochondria") | is.na(Family))
V4.no.mito # 5136 taxa and 176 samples
V4.no.organelles <- subset_taxa(V4.no.mito, (Order != "Chloroplast") | is.na(Order))
V4.no.organelles # 5128 taxa and 176 samples

# Get final sample sums (i.e. sequencing depth per sample) after these steps
sample_data(V4.no.organelles)$sample_sums_final <- sample_sums(V4.no.organelles)



# Add species assignment to ASVs with 100% match to references
# First extract taxonomy table from phyloseq object
taxa <- tax_table(V4.no.organelles)
taxa <- addSpecies(taxa, here("data", "references", "silva_species_assignment_v138.fa.gz"))
# 
# # How many have a species assignment?
length(taxa[, 7]) - sum(is.na(taxa[, 7]))
# # 183 out of 5128 ASVs have a species assignment (3.6%)
# 
saveRDS(taxa, here("data", "process", "taxonomy_with_species_assignment.rds"))
# taxa <- readRDS(here("data", "process", "taxonomy_with_species_assignment.rds"))

# Merging this directly with the phyloseq object doesn't work so need to
# extract the components and make the phyloseq object again
metadata <- sample_data(V4.no.organelles)
seqtab <- otu_table(V4.no.organelles)

# Combine into phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

# Change ASV names ####
# With the combined phyloseq object, change the ASV names from the full ~400 bp
# DNA sequence to much shorter and more manageable names ASVx, where x is the
# ASV number in decreasing order of overall abundance. The ASV DNA sequences 
# are stored in a separate object called `dna`, kept with the phyloseq object.
seqs <- Biostrings::DNAStringSet(taxa_names(ps))
names(seqs) <- taxa_names(ps)
ps <- merge_phyloseq(ps, seqs)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Write pruned and filtered phyloseq object to file 
saveRDS(ps, here("data", "process", "phyloseq_after_pruning_filtering.rds"))
