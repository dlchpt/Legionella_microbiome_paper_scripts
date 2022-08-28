# Figure: Taxonomy strip charts

library(ggpubr)
library(here)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)

# Load plot settings and final data sets
source(here("code", "03_settings_and_data_for_figures.R"))


# Panel A: Phylum level ####
# Agglomerate subsampled RA object to phylum level, keep NAs here
V4.phylum <- V4.RA %>% tax_glom("Phylum", NArm = FALSE)

# Keep only those with a mean abundance > 0.1%
V4.main.phylum <- filter_taxa(V4.phylum, function(x) mean(x) > 0.001, TRUE)

V4.phylum.melted <- as_tibble(phyloseq::psmelt(V4.main.phylum), rownames = "Rows") %>%
  mutate(Abundance_adj = na_if(Abundance, '0')) %>%
  mutate(City = factor(City, 
                       levels = c("Athens", "Copenhagen", "Rome", "Warsaw"))) %>%
  mutate(Phylum = replace_na(Phylum, "unclassified")) %>%
  filter(Phylum != "unclassified") %>%
  mutate(Phylum = fct_reorder(as_factor(Phylum), Abundance, mean, .desc = TRUE))

top.phyla <- levels(V4.phylum.melted$Phylum)[1:10]  

phylum.plot <- V4.phylum.melted %>%
  filter(Phylum %in% top.phyla) %>%
  # mutate(Phylum = fct_relevel(Phylum, "unclassified", after = Inf)) %>%
  ggplot(aes(x = City, y = Abundance_adj * 100, colour = City)) +
  geom_boxplot(colour = "black") +
  geom_jitter(aes(shape = CFU_pres_abs), width = 0.2, alpha = 0.8) +
  facet_wrap(~ Phylum, ncol = 5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) +
  rotate_x_text(angle = 45) +
  scale_colour_manual(values = cbPalette.cities, na.value = "light grey") +
  scale_shape_manual(name = "Legionella cult.", values = shapes.cfu, na.value = 1) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  rremove("legend") +
  rremove("x.text") +
  labs(x = "", y = "Relative abundance (%)")



# Panel B: Genus level ####
# Agglomerate subsampled RA object to genus level, keep NAs here
V4.genus <- V4.RA %>% tax_glom("Genus", NArm = FALSE)

# Keep only those with a mean abundance > 0.1%
V4.main.genus <- filter_taxa(V4.genus, function(x) mean(x) > 0.001, TRUE)

V4.genus.melted <- as_tibble(phyloseq::psmelt(V4.main.genus), rownames = "Rows") %>%
  mutate(Abundance_adj = na_if(Abundance, '0')) %>%
  mutate(City = factor(City, 
                       levels = c("Athens", "Copenhagen", "Rome", "Warsaw"))) %>%
  mutate(Genus = replace_na(Genus, "unclassified")) %>%
  filter(Genus != "unclassified") %>%
  mutate(Genus = fct_reorder(as_factor(Genus), Abundance, mean, .desc = TRUE))

top.genera <- levels(V4.genus.melted$Genus)[1:15]

genus.plot <- V4.genus.melted %>%
  filter(Genus %in% top.genera) %>%
  ggplot(aes(x = City, y = Abundance_adj * 100, colour = City)) +
  geom_boxplot(colour = "black") +
  geom_jitter(aes(shape = CFU_pres_abs), width = 0.2, alpha = 0.8) +
  facet_wrap(~ Genus, ncol = 5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) +
  rotate_x_text(angle = 45) +
  scale_colour_manual(values = cbPalette.cities, na.value = "light grey") +
  scale_shape_manual(name = "Legionella cult.", values = shapes.cfu, na.value = 1) +
  rremove("legend") +
  rremove("x.text") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(x = "", y = "Relative abundance (%)")


# Combined taxonomy figure
tax.fig <- ggarrange(phylum.plot, NULL, genus.plot,
          ncol = 1,
          heights = c(1, 0.1, 1.5),
          labels = c("A", "", "B"))

ggsave(plot = tax.fig, filename = "results/figures/Figure1_taxonomy_v2.pdf", width=9, height=9)




# Supplemental figure: other taxonomic levels ####

# Class level ####
# Agglomerate subsampled RA object to class level, keep NAs here
V4.class <- V4.RA %>% tax_glom("Class", NArm = FALSE)

# Keep only those with a mean abundance > 0.01%
V4.main.class <- filter_taxa(V4.class, function(x) mean(x) > 0.0001, TRUE)

V4.class.melted <- as_tibble(phyloseq::psmelt(V4.main.class), rownames = "Rows") %>%
  mutate(Abundance_adj = na_if(Abundance, '0')) %>%
  mutate(City = factor(City, 
                       levels = c("Athens", "Copenhagen", "Rome", "Warsaw"))) %>%
  mutate(Class = replace_na(Class, "unclassified")) %>%
  mutate(Class = fct_reorder(as_factor(Class), Abundance, mean, .desc = TRUE)) %>%
  filter(Class != "unclassified")

top.classes <- levels(V4.class.melted$Class)[1:15]

class.plot <- V4.class.melted %>%
  filter(Class %in% top.classes) %>%
  ggplot(aes(x = City, y = Abundance_adj * 100, colour = City)) +
  geom_boxplot(colour = "black") +
  geom_jitter(aes(shape = CFU_pres_abs), width = 0.2, alpha = 0.8) +
  facet_wrap(~ Class, ncol = 5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) +
  rotate_x_text(angle = 90) +
  scale_colour_manual(values = cbPalette.cities, na.value = "light grey") +
  scale_shape_manual(name = "Legionella cult.", values = shapes.cfu, na.value = 1) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  rremove("legend") +
  labs(title = "Main bacterial classes", x = "", y = "Relative abundance (%)")

ggsave(plot = class.plot, filename = "results/figures/FigureSX_taxonomy_class.pdf", width=9, height=6)


# Nitrospira only ####
V4.nitro <- subset_taxa(V4.RA, Genus == "Nitrospira")

# Keep only those with a mean abundance > 0.01%
V4.main.nitro <- filter_taxa(V4.nitro, function(x) mean(x) > 0.0001, TRUE)

V4.nitro.melted <- as_tibble(phyloseq::psmelt(V4.main.nitro), rownames = "Rows") %>%
  mutate(Abundance_adj = na_if(Abundance, '0')) %>%
  mutate(City = factor(City, 
                       levels = c("Athens", "Copenhagen", "Rome", "Warsaw"))) %>%
  mutate(OTU = fct_reorder(as_factor(OTU), Abundance, mean, .desc = TRUE))

top.nitro <- levels(V4.nitro.melted$OTU)[1:20]

nitro.plot <- V4.nitro.melted %>%
  filter(OTU %in% top.nitro) %>%
  ggplot(aes(x = City, y = Abundance_adj * 100, colour = City)) +
  geom_boxplot(colour = "black") +
  geom_jitter(aes(shape = CFU_pres_abs), width = 0.2, alpha = 0.8) +
  facet_wrap(~ OTU, ncol = 5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) +
  rotate_x_text(angle = 90) +
  scale_colour_manual(values = cbPalette.cities, na.value = "light grey") +
  scale_shape_manual(name = "Legionella cult.", values = shapes.cfu, na.value = 1) +
  rremove("legend") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(title = "Main Nitrospira ASVs", x = "", y = "Relative abundance (%)")


ggsave(plot = nitro.plot, filename = "results/figures/FigureSX_taxonomy_nitrospira_ASVs.pdf", width=9, height=8)



## Count number of ASVs in each phylum and genus
# Summarised data for editing
phylum.summary <- as_tibble(phyloseq::psmelt(V4.RA), rownames = "Rows") %>%
  select(Phylum, OTU) %>%
  distinct() %>%
  group_by(Phylum) %>%
  summarise(no.asvs = n()) %>%
  filter(Phylum %in% top.phyla)

genus.summary <- as_tibble(phyloseq::psmelt(V4.RA), rownames = "Rows") %>%
  select(Genus, OTU) %>%
  distinct() %>%
  group_by(Genus) %>%
  summarise(no.asvs = n()) %>%
  filter(Genus %in% top.genera)
