library(DEGreport)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(readr)
library(readxl)
library(microbiome)
library(phyloseq)
library(ggplot2)
library(writexl)
library(vegan)
library(here)

# Define age and measure type
which_age <- "21"
measure_type <- "Digesta"

# Read and select metadata
meta <- read_excel("../../data/metadata/metadata.xlsx")
selected_meta <- meta %>%
  filter(TYPE %in% measure_type,
         GUT_SECTION == "Cecum",
         TREATMENT %in% c("PFA", "AGP", "CTR"),
         AGE == which_age)

# Set sample-id as rownames and order by row names
selected_meta <- selected_meta %>%
  column_to_rownames(var = "sample-id") %>%
  arrange(row.names(.))
selected_meta <-select(selected_meta,TREATMENT)

# Read species data 
species_data <- read_excel("../../data/processed/new_GTDB_cecum_species_raw_abundance_filterd_low.xlsx") %>%
  column_to_rownames(var = "sample-id")

# Select and reorder species data based on metadata
selected_species_data <- species_data %>%
  filter(rownames(species_data) %in% rownames(selected_meta)) %>%
  arrange(row.names(.))

# Perform NMDS analysis
m_selected_species_data <- as.matrix(selected_species_data)
nmds <- metaMDS(m_selected_species_data, distance = "bray")

# Fit environmental variables to NMDS
en <- envfit(nmds, selected_meta, permutations = 1000, na.rm = TRUE)

# Extract NMDS scores and environmental factor coordinates
data.scores <- as.data.frame(scores(nmds)$sites) %>%
  mutate(Treatment = selected_meta$TREATMENT)

en_coord_cat <- as.data.frame(scores(en, "factors"))
view(data.scores)

# Plot the NMDS result
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Treatment), size = 3, alpha = 0.7) + 
  scale_colour_manual(values = c("steelblue", "grey", "gold"))  + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Treatment")
ggsave(paste0("../../results/figures/envfit_", measure_type, "_", which_age, "_without_mean.png"), width = 3, height = 2)
dev.off()

