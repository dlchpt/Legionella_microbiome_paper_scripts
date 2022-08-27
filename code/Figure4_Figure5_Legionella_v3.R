# Figure 1: Legionella prevalence and abundance
# version 2: New panel, heatmap of Legionella ASVs grouped by city and building
# version 3: New panel, count of samples per culture category
# Culturable vs 16S data, proportion L.pneumophila, heatmap by city and building

library(ggpubr)
library(ggtext)
library(here)
library(pheatmap)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(UpSetR)
library(viridis)

# Load plot settings and final data sets
source(here("code", "settings_and_data_for_figures.R"))

lod <- 1/nseqs*100

city.levels <- c("Athens", "Copenhagen", "Rome", "Warsaw")
cfu.levels <- c("0", "10<sup>1</sup>", "10<sup>2</sup>", "10<sup>3</sup>",
                "10<sup>4</sup>", "10<sup>5</sup>", "10<sup>6</sup>+")

# Panel A: Count per CFU category ####
fig1a <- samples_V4 %>%
  mutate(City = factor(City, levels = city.levels)) %>%
  filter(!is.na(CFU_cat)) %>%
  ggplot(aes(x = CFU_cat, fill = City)) +
  geom_bar() +
  facet_wrap(~ City, ncol = 4) +
  scale_fill_manual(values = cbPalette.cities) +
  scale_x_discrete(labels = cfu.levels) +
  labs(x = "", 
       y = "Sample count") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
        axis.text.x = element_markdown())
  

# Panel B: Proportion Legionella genus in 16S data vs CFU category ####
fig1b <- samples_V4 %>%
  mutate(City = factor(City, levels = city.levels)) %>%
  filter(!is.na(CFU_cat)) %>%
  ggplot(aes(x = CFU_cat, y = (Leg_genus_sum*100 + 2/3*lod), colour = City, shape = CFU_pres_abs)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(0.01, 0.1, 1, 10, 100)) +
  facet_wrap(~ City, ncol = 4) +
  scale_colour_manual(values = cbPalette.cities) +
  scale_shape_manual(values = shapes.cfu, na.value = 8) +
  geom_hline(yintercept = lod, colour = "black", linetype = "dashed") +
  scale_x_discrete(labels = cfu.levels) +
  labs(x = "Legionella culture counts (CFU/L category)", 
       y = "Legionella % in 16S data ") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_markdown())
  

fig1ab <- ggarrange(fig1a, fig1b,
          nrow = 2,
          labels = c("A", "B"),
          heights = c(1, 1),
          align = "v")


# fig1ab

# Panel C: Proportion of 16S Legionella seqs that are L. pneumo ####
fig1c.plot <- samples_V4 %>%
  mutate(City = factor(City, levels = city.levels)) %>%
  filter(Leg_genus_sum > 0 & !is.na(CFU_pres_abs)) %>%
  mutate(pneumo_prop = Leg_pneumo_sum / Leg_genus_sum) %>%
ggplot(aes(x = City, y = pneumo_prop*100, colour = City, shape = CFU_pres_abs)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(size = 3, position = position_jitterdodge(), alpha = 0.7) +
  scale_colour_manual(values = cbPalette.cities) +
  scale_shape_manual(values = shapes.cfu, na.value = 8) +
  labs(x = "", y = "*L.pneumo* % (of Legionella 16S data)") +
  theme(legend.position = "none",
        axis.title.y = element_markdown())

fig1c <- ggarrange(NULL, fig1c.plot, NULL,
                   ncol = 1,
                   heights = c(0.25, 1, 0.25),
                   labels = c("", "C", ""))


# Figure 4: Heatmap of Legionella ASVs ####

# Set prevalence filter (ASV must occur in more than this number of samples)
prev.threshold <- 4

# Get phyloseq object of just Legionella
V4.leg.ps <- V4.RA %>% subset_taxa(Genus == "Legionella")

# Prevalence filter: keep only Legionella ASVs above a set prevalence threshold
V4.leg.filt <- filter_taxa(V4.leg.ps, function(x) sum(x > 0) > prev.threshold, TRUE)


# Need to calculate intersections to be able to manually set the order of rows in the heatmap
# Melt the prevalence-filtered phyloseq object
V4.leg.genus <- V4.leg.filt %>% psmelt()

# Group by city and OTU
V4.leg.genus.grouped <- V4.leg.genus %>%
  group_by(City, OTU) %>%
  summarise(n = n(),
            Tot.with.OTU = sum(Abundance > 0)) %>%
  filter(Tot.with.OTU > 1) #%>% # OTUs must occur twice in a city to be included in that city's list
  # mutate(OTU = factor(OTU, levels = leg.order))

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

intersect.leg.asvs

# Then get just the ASVs unique to each city, and those found in more than one city.
# These can be computed from the New_data object created during the intersect plot generation.

# Get the list of ASVs in the order they're inputted into intersect command
intersect.asvs.input <- unique(unlist(leg.all.list, use.names = FALSE))

# Get the pres/abs table from intersect object, add ASVs and columns with sums
pres.abs <- intersect.leg.asvs$New_data %>%
  bind_cols(ASV = intersect.asvs.input) %>%
  mutate(total = Athens + Copenhagen + Rome + Warsaw,
         ASV.adj = as.double(gsub("ASV", "", ASV))) %>%
  arrange(desc(total))

# Need to order the ASVs with total > 1 by their ASV number, but keep the order
# of the City-specific ASVs the same...

# First, split the dataframe into city-specific and across multiple cities
# City-specific: in the desired order so no need to change anything
pres.abs.city <- pres.abs %>%
  filter(total == 1)

# Across cities: need to order by total then by ASV.adj
pres.abs.across.cities <- pres.abs %>%
  filter(total > 1) %>%
  arrange(desc(total), ASV.adj)

# Finally, paste these two ordered tables back together and pull out the list of ASVs
# to be used for future ordering
ordered.leg.asvs <- bind_rows(pres.abs.across.cities, pres.abs.city) %>%
  pull(ASV)


# Now that I have a specified order for the ASVs, can proceed with heatmap!
# Extract OTU table from Legionella phyloseq object, with ASVs as cols and samples as rows
# log2 transform, then convert to tibble
otus.leg.raw <- otu_table(V4.leg.filt)
otus.leg.log <- log2(otus.leg.raw + 0.000001)

# Replace log2(0.000001) with NA, since 0.00001 was just added to zeros to permit log transform
otus.leg <- na_if(otus.leg.log, log2(0.000001))

otus.leg.tib <- as_tibble(as.data.frame(otus.leg), rownames = "Sample")


# Extract metadata from Legionella-only phyloseq object
metadata.leg <- as_tibble(sample_data(V4.leg.filt), rownames = "Sample") %>%
  select(Sample, City, Bldg_ID, Location_in_building, 
         CFU_L, CFU_cat, CFU_pres_abs, Temp_degC) %>%
  mutate(Bldg_no = str_sub(Bldg_ID, -2, -1))

metadata.leg$City <- fct_relevel(metadata.leg$City, 
                                 levels = c("Athens", "Copenhagen", "Rome", "Warsaw"))

# Join log-transformed OTU table and metadata for arranging samples into desired order
otus.leg.with.metadata <- left_join(otus.leg.tib, metadata.leg) %>%
  arrange(City, Bldg_no, Location_in_building)

# Now that the samples are listed in the desired order, extract the ASV table
# and metadata again, and transpose the ASV table
otus.leg.sorted <- otus.leg.with.metadata %>%
  select(starts_with("ASV")) %>%
  as.matrix()

rownames(otus.leg.sorted) <- otus.leg.with.metadata$Sample
otus.leg.sorted.t <- t(otus.leg.sorted)


# Pull out named species associated with specific ASVs, rename ASV280
leg.sp <- V4.leg.genus %>%
  select(OTU, Species) %>%
  distinct() %>%
  mutate(Species = if_else(OTU == "ASV280", "rubrilucens.taurinensis", Species)) %>%
  rename(ASV = OTU) %>%
  filter(!is.na(Species))

# The current transposed matrix has ASV names as rownames. Need to sort these
# based on the previously generated character vector ordered.leg.asvs
otus.trans.tib <- as_tibble(otus.leg.sorted.t, rownames = "ASV")

otus.trans.tib.ord <- left_join(tibble("ASV" = ordered.leg.asvs), otus.trans.tib) %>%
  left_join(leg.sp) %>%
  unite(col = "ASV.adj", ASV, Species, sep = "_", remove = FALSE) %>%
  mutate(ASV.adj = gsub("_NA", "", ASV.adj))

otus.trans.ord <- otus.trans.tib.ord %>%
  select(-ASV, -ASV.adj, -Species) %>%
  as.matrix()

rownames(otus.trans.ord) <- otus.trans.tib.ord$ASV.adj

# Prepare column annotation dataframe from photo metadata tibble
column.annotations <- data.frame(City = metadata.leg$City,
                                 Bldg = metadata.leg$Bldg_no,
                                 Leg.pres.abs = metadata.leg$CFU_pres_abs,
                                 row.names = metadata.leg$Sample)

# Prepare annotation colours
viridis.3 <- viridis(3, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
viridis.bldgs <- c("B1" = viridis.3[1], "B2" = viridis.3[2], "B3" = viridis.3[3])

annotation.colours <- list(City = cbPalette.cities,
                           Bldg = viridis.bldgs,
                           Leg.pres.abs = cbPalette.CFU.pres.abs)


# Manually set the column gaps
City_col_gaps <- c(43, 86, 130)

ASV_row_gaps <- c(10, 18, 31, 38)


# Finally, draw the heatmap!
# Rows clustered by Euclidean distance
fig1d <- pheatmap(
  mat               = otus.trans.ord,
  color             = viridis(10),
  border_color      = NA,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = column.annotations,
  annotation_colors = annotation.colours,
  drop_levels       = TRUE,
  gaps_col          = City_col_gaps,
  gaps_row          = ASV_row_gaps,
  na_col            = "grey"
)

fig1d
# Only ASVs present in at least 5 samples overall are shown.
# An ASV is included on a City list if it occurs in at least two samples from that city.

# Darn it, pheatmap objects cannot easily be included in ggarrange...
# Might cut my losses and add this panel in Illustrator.

ggsave(plot = fig1d, filename = "results/figures/Figure1c_Legionella_heatmap.pdf", width=12, height=7)



# Combined figure, panels A and B
fig1 <- ggarrange(fig1ab, fig1c,
                  nrow = 1,
                  labels = c("", ""),
                  widths = c(2, 1))

ggsave(plot = fig1, filename = "results/figures/Figure1_Legionella.pdf", width=12, height=6)
