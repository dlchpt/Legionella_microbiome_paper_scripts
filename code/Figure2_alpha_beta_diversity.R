# Figure 2: Alpha and beta diversity
# Population richness, shannon, PCA

library(breakaway)
library(cowplot)
library(factoextra)
library(ggpubr)
library(grid)
library(here)
library(microbiome)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(vegan)
library(zCompositions)

# Load plot settings and final data sets
source(here("code", "03_settings_and_data_for_figures.R"))

## Multipanel figure: Richness, shannon, PCA ####
# Panel A: Breakaway richness estimate ####
# Breakaway calculation - takes ages so will comment out once run
br_V4 <- breakaway(V4)

# Save breakaway output to load quickly without having to re-run
saveRDS(br_V4, here("data", "process", "16S_breakaway.rds"))

br_V4 <- readRDS(here("data", "process", "16S_breakaway.rds"))

# Gather results into tibble with metadata
br_V4_ord <- br_V4 %>% 
  summary %>%
  full_join(samples_V4, by = c("sample_names" = "SampleID")) %>%
  mutate(City = fct_relevel(City, c("Athens", "Copenhagen", "Rome", "Warsaw"))) %>%
  dplyr::arrange(City, Season_year, Location_in_building) %>%
  mutate(sample_names = fct_inorder(sample_names)) # set the order of rows

rm(br_V4)

fig2a <- br_V4_ord %>%
  filter(!is.na(CFU_pres_abs)) %>%
  ggplot(aes(x = City, y = estimate, colour = City, shape = CFU_pres_abs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge(), alpha = 0.7) +
  scale_colour_manual(values = cbPalette.cities) +
  scale_shape_manual(values = shapes.cfu) +
  theme(legend.position = "none") +
  rotate_x_text() +
  labs(x = "", y = "Est. pop. richness")


# Panel B: Shannon index ####
# For this, I use the rarefied data set (but not relative abundance)
# microbiome package function to compute the shannon index
diver <- microbiome::diversity(V4.rare, index = "shannon") %>%
  as_tibble(rownames = "SampleID") %>%
  left_join(samples_V4, by = "SampleID") %>%
  mutate(City = factor(City, levels = c("Athens", "Copenhagen", "Rome", "Warsaw")))


fig2b <- diver %>%
  filter(!is.na(CFU_pres_abs)) %>%
  ggplot(aes(x = City, y = shannon, colour = City, shape = CFU_pres_abs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge(), alpha = 0.7) +
  scale_colour_manual(values = cbPalette.cities) +
  scale_shape_manual(values = shapes.cfu) +
  theme(legend.position = "none") +
  rotate_x_text() +
  labs(x = "", y = "Shannon index")



# Panel C: PCA on Aitchison distance ####
# Remove low-abundance reads, which don't contribute to the PCA
V4.trim.unique.main <- prune_taxa(taxa_sums(V4) > 19, V4)

otus_V4 <- as.data.frame(t(otu_table(V4.trim.unique.main))) %>%
  as_tibble(rownames = "OTU_ID")

# Convert otu table to matrix and use ASV names as row names
V4_m <- as.matrix(otus_V4[,-1])
rownames(V4_m) <- otus_V4$OTU_ID

# Zero count replacement
# this function expects the samples to be in rows and OTUs to be in columns
# so the matrix is transposed on input, and then back again on output
V4_mz <- t(zCompositions::cmultRepl(t(V4_m), method="CZM", output="p-counts"))

# convert to proportions by sample (columns) using the apply function
V4_prop <- apply(V4_mz, 2, function(x){x/sum(x)})

# make our compositional dataset with the center mean log ratio function
# transpose the output in preparation for the PCA, which needs samples
# in rows
V4_comp <- t(apply(V4_prop, 2, function(x){log(x) - mean(log(x))}))

# Run the PCA on the CLR-transformed data. Don't need to scale as the CLR
# transform has already done that
V4_PCA <- prcomp(V4_comp, scale = FALSE)

# Make a vector for the CFU_pres_abs symbol
samples_V4_pca <- samples_V4 %>%
  mutate(CFU_shape = case_when(
    CFU_pres_abs == "neg" ~ shapes.cfu["neg"],
    CFU_pres_abs == "pos" ~ shapes.cfu["pos"],
    is.na(CFU_pres_abs) ~ 1)) %>%
  mutate(CFU_pres_abs = replace_na(CFU_pres_abs, "not.measured"))

# Figure: main ordination
# manual.xlim <- c(-50, 100)
# manual.ylim <- c(-110, 50)
manual.xlim <- c(-125, 125)
manual.ylim <- c(-125, 75)

pmain_V4 <- fviz_pca_biplot(V4_PCA, 
                            axes = c(1, 2), # which PCs to plot
                            repel = TRUE, # avoid text overlapping
                            # Samples (i.e. individuals)
                            geom = "point", 
                            col.ind = "white",
                            pointsize = 0,
                            # pointshape = samples_V4_pca$CFU_shape, alpha = 0.7,
                            mean.point = FALSE,
                            # ASVs (i.e. variables)
                            geom.var = c("arrow", "text"),
                            select.var = list(contrib=12),
                            col.var = "black") + 
  # coord_fixed(1) + # scale the coordinates based on % variance explained
  coord_cartesian(xlim = manual.xlim, ylim = manual.ylim) +
  geom_point(aes(colour = samples_V4_pca$City, shape = as.factor(samples_V4_pca$CFU_pres_abs))) +
  scale_colour_manual(name = "City", values = cbPalette.cities) +
  scale_shape_manual(values = c("not.measured" = 1, "neg" = 8, "pos" = 19)) +
  labs(title = "", shape = "Cultured Legionella")

# pmain_V4

# Marginal boxplots
# Extract results of PCA as lists
V4_PCA_inds <- get_pca_ind(V4_PCA)

# Make a data frame from just the coordinates for PC1/2/3
V4_PCA_coords <- as_tibble(V4_PCA_inds$coord[,1:3], rownames = "SampleID")

# Join with streamlined metadata
V4_PCA_coords_metadata <- left_join(V4_PCA_coords, samples_V4_pca, by = "SampleID")


# Marginal boxplot of x (top panel) and y (right panel)
xplot_V4 <- ggboxplot(V4_PCA_coords_metadata, x = "City", y = "Dim.1", 
                      color = "City", fill = "City",
                      alpha = 0.5, ggtheme = theme_bw()) +
  scale_colour_manual(name = "City", values = cbPalette.cities) +
  scale_fill_manual(name = "City", values = cbPalette.cities) +
  # coord_cartesian(ylim = c(-40, 60)) +
  coord_flip(ylim = manual.xlim) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


yplot_V4 <- ggboxplot(V4_PCA_coords_metadata, x = "City", y = "Dim.2",
                      color = "City", fill = "City",
                      alpha = 0.5, ggtheme = theme_bw()) +
  scale_colour_manual(name = "City", values = cbPalette.cities) +
  scale_fill_manual(name = "City", values = cbPalette.cities) +
  coord_cartesian(ylim = manual.ylim) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Cleaning the plots
pmain_V4 <- pmain_V4 + rremove("legend")
yplot_V4 <- yplot_V4 + clean_theme() + rremove("legend")
xplot_V4 <- xplot_V4 + clean_theme() + rremove("legend")


p1_box_V4 <- insert_xaxis_grob(pmain_V4, xplot_V4, position = "top")
p2_box_V4 <- insert_yaxis_grob(p1_box_V4, yplot_V4, position = "right")
fig2c <- ggdraw(p2_box_V4)
# fig2c



# Combined figure ####
row1 <- ggarrange(fig2a, NULL, fig2b, 
                                     ncol = 3, 
                                     labels = c("A", "", "B"), 
                                     widths = c(1, 0.1, 1),
                                     align = "h")

row2 <- ggarrange(fig2c, ncol = 1, labels = c("C"))

fig2 <- ggarrange(row1, NULL, row2,
                     nrow = 3,
                     heights = c(1, 0.01, 1.8))

# fig2

ggsave(plot = fig2, filename = "results/figures/Figure2_alpha_beta_diversity.pdf", width=9, height=9)
