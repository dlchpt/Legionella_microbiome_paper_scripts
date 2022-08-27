# Figure 3: Intersection plot of Legionella ASVs and all ASVs

library(ggpubr)
library(here)
library(microbiome)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(UpSetR)

# Load plot settings and final data sets
source(here("code", "settings_and_data_for_figures.R"))


# Panel A: Intersection of all Legionella ASVs ####
# Select only ASVs classified as Legionella
V4.leg.genus <- V4.RA %>% subset_taxa(Genus == "Legionella") %>%
  psmelt()

# Set prevalence filter (ASV must occur in more than this number of samples)
prev.threshold <- 0 # want to keep all ASVs here


V4.leg.genus.grouped <- V4.leg.genus %>%
  group_by(City, OTU) %>%
  summarise(n = n(),
            Tot.with.OTU = sum(Abundance > 0),
            Prevalence = sum(Abundance > 0) / n(),
            Mean.abund = mean(Abundance),
            Tot.abund = sum(Abundance)) %>%
  filter(Tot.with.OTU > prev.threshold)

# Pull out a character vector for each city, with the names of Legionella ASVs detected there
leg.athens <- V4.leg.genus.grouped %>%
  filter(City == "Athens") %>% pull(OTU) %>% as.character()

leg.rome <- V4.leg.genus.grouped %>%
  filter(City == "Rome") %>% pull(OTU) %>% as.character()

leg.warsaw <- V4.leg.genus.grouped %>%
  filter(City == "Warsaw") %>% pull(OTU) %>% as.character()

leg.copenh <- V4.leg.genus.grouped %>%
  filter(City == "Copenhagen") %>% pull(OTU) %>% as.character()

# Make a list of all four city Legionella vectors
leg.all.list <- list(Athens = leg.athens,
                     Copenhagen = leg.copenh,
                     Rome = leg.rome,
                     Warsaw = leg.warsaw)

# Make a list of the intersections to highlight in colour
# Common core in orange, region-specific in blue
queries.list <- list(
  list(query = intersects, params = list("Athens", "Copenhagen", "Rome", "Warsaw"), color = "black", active = T),
  list(query = intersects, params = list("Athens"), color = cbPalette.cities["Athens"], active = T),
  list(query = intersects, params = list("Copenhagen"), color = cbPalette.cities["Copenhagen"], active = T),
  list(query = intersects, params = list("Rome"), color = cbPalette.cities["Rome"], active = T),
  list(query = intersects, params = list("Warsaw"), color = cbPalette.cities["Warsaw"], active = T)
)

# Make figure and save as pdf
intersect.leg.asvs <- upset(fromList(leg.all.list), main.bar.color = "dark grey", order.by = "freq", sets.x.label = "Total Leg. ASVs", queries = queries.list)

# intersect.leg.asvs

pdf(file=here("results", "figures", "Figure3a_Legionella_ASV_intersection.pdf"), width = 5, height = 4)
intersect.leg.asvs
dev.off()





# Panel B: Intersection of all ASVs ####
# First, get the list of ASVs from each city. For this, I'll use the phyloseq
# object that is subsampled to equal size - otherwise the sequencing batch
# effect could skew this analysis

# Vector of region names
regions <- unique(samples_V4$City)

for (i in regions){
  # Subsample to each region
  V4.region <- subset_samples(V4.rare, City == i)
  V4.region <- prune_taxa(taxa_sums(V4.region) > 0, V4.region)
  all.list <- taxa(V4.region)
  assign(paste0("all.list.V4.", i), all.list)
}

# Make a list of the named vectors, one per city
all.list.V4.regions <- list(Athens = all.list.V4.Athens,
                            Copenhagen = all.list.V4.Copenhagen,
                            Rome = all.list.V4.Rome,
                            Warsaw = all.list.V4.Warsaw)

# Make a list of the intersections to highlight in colour
# Common core in orange, region-specific in blue
queries.list <- list(
  list(query = intersects, params = list("Athens", "Copenhagen", "Rome", "Warsaw"), color = "black", active = T),
  list(query = intersects, params = list("Athens"), color = cbPalette.cities["Athens"], active = T),
  list(query = intersects, params = list("Copenhagen"), color = cbPalette.cities["Copenhagen"], active = T),
  list(query = intersects, params = list("Rome"), color = cbPalette.cities["Rome"], active = T),
  list(query = intersects, params = list("Warsaw"), color = cbPalette.cities["Warsaw"], active = T)
)

# Make figure and save as pdf
intersect.asvs <- upset(fromList(all.list.V4.regions), main.bar.color = "dark grey", order.by = "freq", empty.intersections = "on", sets.x.label = "Total richness (no. ASVs)", queries = queries.list)

pdf(file=here("results", "figures", "Figure3b_ASV_intersection.pdf"), width = 6, height = 4)
intersect.asvs
dev.off()

