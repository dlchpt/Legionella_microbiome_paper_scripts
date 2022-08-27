# FIGURES: Settings and datasets
# This is maintain consistency across figures, and to make it easy to change
# them all without having to adapt multiple different scripts.

# I'll also load the final pruned and filtered dataset, recode the CFU/L data
# and other variables, and compute the rarefied, relative abundance, and
# taxonomy glommed phyloseq objects, since these are used in multiple figures.

# Required libraries ####
library(here); packageVersion("here")
library(microbiome); packageVersion("microbiome")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")

theme_set(theme_classic())
set.seed(65536)


# Plot settings ####

## Colour palettes
cbPalette.cities <- c(
  "Warsaw" = "#D55E00", # red
  "Copenhagen" = "#E69F00", # orange
  "Rome" = "#0072B2", # dark blue
  "Athens" = "#009E73" # green
  )


cbPalette.CFU.pres.abs <- c("neg" = "grey",
                            "pos" = "dark red")

cbPalette.CFU.cat <- c("#5E4FA2", "#48A0B2", "#A1D9A4", "#EDF7A3", "#FEE899", 
                       "#FBA45C", "#E25249", "#9E0142")

cbPalette.class <- c(
  "Gammaproteobacteria" = "#D73027", # red Proteobact
  "Alphaproteobacteria" = "#FDAE61", # mid orange Proteobact
  "Bacilli" = "#A6D96A", # light green Firmicutes
  "Deinococci" = "#1A9850", # dark green Deino
  "Nitrospiria" = "#F0E442", # yellow Nitrospirota
  "Bacteroidia" = "#762A83", # dark purple Bacteroi
  "Blastocatellia" = "#74ADD1", # mid blue Acidobact
  "Acidobacteriae" = "#ABD9E9", # light blue Acidobact
  "Actinobacteria" = "#C51B7D", # dark pink Actino
  "Thermoleophilia" = "#F1B6DA", # light pink Actino
  "Acidimicrobiia" = "#C2A5CF", # light purple Actino
  "Gemmatimonadetes" = "#878787", # grey Gemma
  "Polyangia" = "#66BD63", # mid green
  "Bdellovibrionia" = "#000000", # black
  "Phycisphaerae" = "#4575B4", # dark blue
  "Chlamydiae" = "#9970AB" # mid purple
)

## Shapes
shapes.compartments <- c("MW" = 19, "BB" = 17, "HRW" = 15,
                         "FPB" = 8, "AHE" = 13)

shapes.cfu <- c("pos" = 19, "neg" = 8)


# Data for final figures and tables ####
# Full phyloseq object from pruning/filtering script
V4 <- readRDS(here("data", "process", "phyloseq_after_pruning_filtering.rds"))

# To append Legionella rel. abundance, convert the whole object to RA
V4.RA.full <- transform(V4, "compositional")

# Select only ASVs classified as Legionella, then add their total RA to the 
# metadata
V4.leg.genus <- V4.RA.full %>% subset_taxa(Genus == "Legionella")
sample_data(V4)$Leg_genus_sum <- sample_sums(V4.leg.genus)


# Select only ASVs classified as Legionella pneumophila, then add their total 
# RA to the metadata
V4.leg.pneumo <- V4.RA.full %>% subset_taxa(Genus == "Legionella" & 
                                         Species == "pneumophila")
sample_data(V4)$Leg_pneumo_sum <- sample_sums(V4.leg.pneumo)

# Remove temporary phyloseq objects
rm(V4.RA.full, V4.leg.genus, V4.leg.pneumo)

# Extract metadata and recode CFU/L into categorical variables
samples_V4 <- as_tibble(sample_data(V4), rownames = "SampleID") %>%
  mutate(Season_year = factor(Season_year, 
                              levels = c("Winter2017", "Summer2017", 
                                         "Autumn2017", "Winter2018", 
                                         "Spring2018"))) %>%
  mutate(City = as_factor(City)) %>%
  mutate(CFU_cat = case_when(
    CFU_L == 0.0 ~ "none",
    CFU_L > 0.0 & CFU_L <= 10 ~ "10^1",
    CFU_L > 10 & CFU_L <= 100 ~ "10^2",
    CFU_L > 100 & CFU_L <= 1000 ~ "10^3",
    CFU_L > 1000 & CFU_L <= 10000 ~ "10^4",
    CFU_L > 10000 & CFU_L <= 100000 ~ "10^5",
    CFU_L > 10^5 ~ ">10^5",
    TRUE ~ NA_character_)) %>%
  mutate(CFU_pres_abs = case_when(
    CFU_L == 0.0 ~ "neg",
    CFU_L > 0.0 ~ "pos",
    TRUE ~ NA_character_)) %>%
  mutate(CFU_cat = factor(CFU_cat, 
                          levels = c("none", "10^1", "10^2", "10^3", "10^4", 
                                     "10^5", ">10^5", "NA")))


# Make new full phyloseq object with recoded metadata
V4_otus <- otu_table(V4)
V4_tax <- tax_table(V4)
V4_seqs <- refseq(V4)

V4 <- phyloseq(otu_table(V4_otus),
                    sample_data(data.frame(samples_V4, row.names = "SampleID")),
                    tax_table(V4_tax),
                    refseq(V4_seqs))

rm(V4_otus, V4_tax, V4_seqs)

# Modified data sets ####
# Rarefied data for taxonomy plots and plug-in shannon index
nseqs <- 10148
V4.rare <- rarefy_even_depth(V4, sample.size = nseqs, rngseed = 65536)

# Relative abundance of rarefied data
V4.RA <- transform(V4.rare, "compositional")

# Glommed by taxonomy
# V4.class <- V4.RA %>% tax_glom("Class", NArm = FALSE)
# V4.order <- V4.RA %>% tax_glom("Order", NArm = FALSE)
# V4.family <- V4.RA %>% tax_glom("Family", NArm = FALSE)
# V4.genus <- V4.RA %>% tax_glom("Genus", NArm = FALSE)

