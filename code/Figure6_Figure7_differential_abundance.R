# Figure 5: Differential abundance of genera by Legionella CFU

library(ALDEx2)
library(ggpubr)
library(here)
library(microbiome)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)


# Load plot settings and final data sets
source(here("code", "settings_and_data_for_figures.R"))

# Figure 5: ALDEx2 plot at genus level
# This is run from the full phyloseq object, as I need the genus-glommed data
# on raw counts, not on RA
# Exclude samples with no CFU data, keep low-abundance ASVs prior to agglomeration
V4.cfu.genus <- V4 %>% subset_samples(!is.na(CFU_pres_abs)) %>%
  tax_glom("Genus", NArm = FALSE)

# Before removing low-abundance genera, get RA for plot bubble sizes
V4.cfu.genus.RA <- transform(V4.cfu.genus, "compositional")

V4.cfu.genus <- prune_taxa(taxa_sums(V4.cfu.genus) > 49, V4.cfu.genus)

# This takes ages so will save the object and comment out the command once it's run
# aldex2_cfu_genus <- aldex(data.frame(t(otu_table(V4.cfu.genus))), as.character(sample_data(V4.cfu.genus)$CFU_pres_abs), test="t", effect = TRUE, denom="all")

# saveRDS(aldex2_cfu_genus, here("data", "process", "aldex_cfu_genus_pres_abs.rds"))
aldex2_cfu_genus <- readRDS(here("data", "process", "aldex_cfu_genus_pres_abs.rds"))


# Add taxonomic labels to aldex2 output
taxa_info_v4_genus <- data.frame(tax_table(V4.cfu.genus)) %>% 
  rownames_to_column(var = "OTU")

# Trimmed table for display
# Add taxa sums and means, replace Genus NAs with closest classification
sig_aldex2_cfu_genus <- aldex2_cfu_genus %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, effect, wi.eBH) %>%
  left_join(taxa_info_v4_genus) %>%
  left_join(tibble(OTU = names(taxa_sums(V4.cfu.genus.RA)),
                   taxa_sums_RA = taxa_sums(V4.cfu.genus.RA),
                   mean_taxa_RA = taxa_sums_RA / nsamples(V4.cfu.genus.RA))) %>%
  mutate(Genus_adj = case_when(
    is.na(Class) ~ paste0("uncl. phylum ", Phylum),
    is.na(Order) ~ paste0("uncl. class ", Class),
    is.na(Family) ~ paste0("uncl. order ", Order),
    is.na(Genus) ~ paste0("uncl. family ", Family),
    TRUE ~ Genus )) %>%
  mutate(Genus = str_replace(Genus, "Craurococcus-Caldovatus", "Crauro.-Caldo.")) %>%
  mutate(Genus = str_replace(Genus, "CL500-29_marine_group", "CL500-29 marine gr.")) %>%
  filter(!is.na(Genus)) %>%
  dplyr::arrange(Phylum, Class, Order) %>%
  mutate(Class = fct_inorder(Class)) %>%
  mutate(Genus_adj = fct_inorder(Genus_adj)) %>%
  mutate(Genus = fct_inorder(Genus))
  

# Fig 5b ALDEx2 genus ####
fig5b <- sig_aldex2_cfu_genus %>%
  ggplot(aes(x = Genus, y = effect, color = Class)) +
  geom_segment(aes(y = 0, x = Genus, 
                   yend = effect, xend = Genus), 
               color = "grey") +
  geom_point(aes(size = mean_taxa_RA*100), alpha = 1) +
  coord_flip() +
  scale_colour_manual(values = cbPalette.class, na.value = "light grey") +
  labs(y = "Leg. neg <- 0 -> Leg. pos (effect)", 
       x = "", color = "Class", 
       size = "Mean perc.") +
  theme(legend.position = "right",
        panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
  guides(color = guide_legend(reverse = TRUE))




## ALDEx2 at the ASV level
# Exclude samples with no CFU data
V4.cfu <- V4 %>% subset_samples(!is.na(CFU_pres_abs))
V4.cfu <- prune_taxa(taxa_sums(V4.cfu) > 49, V4.cfu)


# This takes ages so will save the object and comment out the command once it's run
# aldex2_cfu <- aldex(data.frame(t(otu_table(V4.cfu))), as.character(sample_data(V4.cfu)$CFU_pres_abs), test="t", effect = TRUE, denom="all")

# saveRDS(aldex2_cfu, here("data", "process", "aldex_cfu_pres_abs.rds"))
aldex2_cfu <- readRDS(here("data", "process", "aldex_cfu_pres_abs.rds"))

# Add taxonomic labels to aldex2 output
# V4
taxa_info_v4 <- data.frame(tax_table(V4.cfu))
taxa_info_v4 <- taxa_info_v4 %>% rownames_to_column(var = "OTU")

# Add taxa sums and means, replace Genus NAs with closest classification
sig_aldex2_cfu <- aldex2_cfu %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, effect, wi.eBH) %>%
  left_join(taxa_info_v4) %>%
  left_join(tibble(OTU = names(taxa_sums(V4.RA)),
                   taxa_sums_RA = taxa_sums(V4.RA),
                   mean_taxa_RA = taxa_sums_RA / nsamples(V4.RA))) %>%
  mutate(Genus_adj = case_when(
    is.na(Class) ~ paste0("uncl. phylum ", Phylum),
    is.na(Order) ~ paste0("uncl. class ", Class),
    is.na(Family) ~ paste0("uncl. order ", Order),
    is.na(Genus) ~ paste0("uncl. family ", Family),
    TRUE ~ Genus )) %>%
  mutate(Genus = str_replace(Genus, "Craurococcus-Caldovatus", "Crauro.-Caldo.")) %>%
  unite("Genus_OTU", OTU, Genus, sep = ":", remove = FALSE, na.rm = FALSE) %>%
  filter(!is.na(Genus)) %>%
  dplyr::arrange(Phylum, Class, Order) %>%
  mutate(Class = fct_inorder(Class)) %>%
  mutate(Genus_OTU = fct_inorder(Genus_OTU)) %>%
  mutate(Genus = fct_inorder(Genus))

# Fig 5a ALDEx2 ASV level ####
fig5a <- sig_aldex2_cfu %>%
  ggplot(aes(x = Genus_OTU, y = effect, color = Class)) +
  geom_segment(aes(y = 0, x = Genus_OTU, 
                   yend = effect, xend = Genus_OTU), 
               color = "grey") +
  geom_point(aes(size = mean_taxa_RA*100), alpha = 1) +
  coord_flip() +
  scale_colour_manual(values = cbPalette.class, na.value = "light grey") +
  labs(y = "Leg. neg <- 0 -> Leg. pos (effect)", 
       x = "", color = "Class", 
       size = "Mean perc.") +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
  guides(color = guide_legend(reverse = TRUE))

# fig5a

fig5 <- ggarrange(fig5a, NULL, fig5b,
                  ncol = 3,
                  widths = c(0.95, 0.1, 1.25),
                  labels = c("A", "", "B"))

# fig5

# Figure 6: Example individual plots ####
# Genus-level plots
# Extract significant genera from genus-agglomerated object
aldex.list.genus <- sig_aldex2_cfu_genus$OTU

V4.aldex.genus <- prune_taxa(aldex.list.genus, V4.cfu.genus.RA)

# Melt the phyloseq object for taxonomic plots
V4.aldex.genus.melted <- as_tibble(phyloseq::psmelt(V4.aldex.genus), rownames = "Rows") %>%
  mutate(Abundance_adj = na_if(Abundance, '0'))


# Same thing at the ASV level
aldex.list <- sig_aldex2_cfu$OTU

V4.aldex <- prune_taxa(aldex.list, V4.RA)
V4.aldex.melted <- as_tibble(phyloseq::psmelt(V4.aldex), rownames = "Rows") %>%
  mutate(Abundance_adj = na_if(Abundance, '0'))



# Limit of detection for figures
lod <- 1/nseqs*100

# Positively associated taxa (all from the individual ASV level, not genus glommed)
sig.pos <- c("Legionella", "Nitrospira", "Craurococcus-Caldovatus", 
             "oc32", "Reyranella")

sig.pos.genus <- c("Bdellovibrio")

pos.genus.table <- V4.aldex.genus.melted %>% 
  filter(Genus %in% sig.pos.genus & !is.na(CFU_pres_abs)) %>%
  dplyr::select(Genus, CFU_pres_abs, Abundance, City, Genus)

pos.asv.table <- V4.aldex.melted %>%
  filter(Genus %in% sig.pos & !is.na(CFU_pres_abs)) %>%
  dplyr::select(Genus, CFU_pres_abs, Abundance, City, Genus)

pos.taxa.fig <- bind_rows(pos.genus.table, pos.asv.table) %>%
  mutate(Genus = recode(Genus, "Legionella" = "Legionella pneumophila ASV142",
                        "Nitrospira" = "Nitrospira ASV23/119/190/249",
                        "oc32" = "Nitrosomonad. oc32 ASV48",
                        "Craurococcus-Caldovatus" = "Craurococcus-Caldovatus ASV31/484",
                        "Reyranella" = "Reyranella ASV91/108",
                        "Bdellovibrio" = "Bdellovibrio (genus)")) %>%
  mutate(Genus = factor(Genus, 
                        levels = c("Legionella pneumophila ASV142", "Nitrospira ASV23/119/190/249",
                                   "Nitrosomonad. oc32 ASV48", "Craurococcus-Caldovatus ASV31/484",
                                   "Reyranella ASV91/108", "Bdellovibrio (genus)"))) %>%
  ggplot(aes(x = CFU_pres_abs, y = (Abundance*100 + 2/3*lod), colour = City)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, position = position_jitterdodge()) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(0.01, 0.1, 1, 10, 100)) +
  scale_colour_manual(values = cbPalette.cities) +
  facet_wrap(~ Genus, nrow = 2) +
  geom_hline(yintercept = lod, colour = "black", linetype = "dashed") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(x = "Culturable Legionella status (pres/abs)", 
       y = "Rel. abund. (%)")

# pos.taxa.fig

# Negatively associated taxa (Staph from ASV level, others from genus glommed)
sig.neg <- "Staphylococcus"
sig.neg.genus <- c("Pseudomonas", "Acinetobacter")
# "Hyphomicrobium", 

neg.genus.table <- V4.aldex.genus.melted %>% 
  filter(Genus %in% sig.neg.genus & !is.na(CFU_pres_abs)) %>%
  dplyr::select(Genus, CFU_pres_abs, Abundance, City, Genus)

neg.asv.table <- V4.aldex.melted %>%
  filter(Genus %in% sig.neg & !is.na(CFU_pres_abs)) %>%
  dplyr::select(Genus, CFU_pres_abs, Abundance, City, Genus)

neg.taxa.fig <- bind_rows(neg.genus.table, neg.asv.table) %>%
  mutate(Genus = recode(Genus, "Staphylococcus" = "Staphylococcus ASV126",
                        "Pseudomonas" = "Pseudomonas (genus)",
                        "Acinetobacter" = "Acinetobacter (genus)")) %>%
  mutate(Genus = factor(Genus, 
                        levels = c("Staphylococcus ASV126", "Pseudomonas (genus)", "Acinetobacter (genus)"))) %>%
  ggplot(aes(x = CFU_pres_abs, y = (Abundance*100 + 2/3*lod), colour = City)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1, position = position_jitterdodge()) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(0.01, 0.1, 1, 10)) +
  scale_colour_manual(values = cbPalette.cities) +
  facet_wrap(~ Genus, ncol = 3) +
  geom_hline(yintercept = lod, colour = "black", linetype = "dashed") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(x = "Culturable Legionella status (pres/abs)", 
       y = "Rel. abund. (%)")

# neg.taxa.fig



# Separate figures! Combined is too busy...

fig6 <- ggarrange(neg.taxa.fig, NULL, pos.taxa.fig,
                        nrow = 3,
                        heights = c(1, 0.1, 2),
                        labels = c("A", "", "B"),
                        common.legend = TRUE,
                        legend = "right")


fig6

ggsave(plot = fig5, filename = "results/figures/Figure5_diff_abund_overall.pdf", width=9, height=9)
ggsave(plot = fig6, filename = "results/figures/Figure6_diff_abund_examples.pdf", width=9, height=6)
