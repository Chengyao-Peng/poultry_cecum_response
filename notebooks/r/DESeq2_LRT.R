## Installing DEGreport on BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

devtools::install_git("https://git@git.bioconductor.org/packages/DEGreport")
# BiocManager::install("DEGreport")
BiocManager::install("apeglm", force=TRUE)

## loading the libraries
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(readr)
library(readxl)
library(microbiome)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(tibble)
library(writexl)

# set the working directory to where the source code is
wd <- getwd()
setwd(wd)

# Define age, measure type and cutoff
padj.cutoff <- 0.05
age_mapping <- c('1' = '3', '2' = '14', '3' = '21', '4' = '35')

# Read and select metadata
meta <- read_excel("../../data/metadata/metadata.xlsx")
selected_meta <- meta %>%
  filter(TYPE == measure_type,
         AGE == day)

# Set sample-id as rownames and order by row names
selected_meta <- selected_meta %>%
  column_to_rownames(var = "sample-id") %>%
  arrange(row.names(.))
selected_meta$TREATMENT <- paste0(selected_meta$TREATMENT)
selected_meta$TREATMENT <- relevel(factor(selected_meta$TREATMENT), ref = "CTR")

# Read count table
data <- read_tsv("../../data/processed/rna_seq_filtered_low.tsv")
data <- data %>% column_to_rownames(var="sample_id")

# Match the count data and the metadata
selected_data <- data[rownames(data) %in% rownames(selected_meta), ]
selected_data <- selected_data[match(rownames(selected_meta), rownames(selected_data)), ]

# Convert the format of the count data
selected_data_T <- t(selected_data)
countDataMatrix <- as.matrix(selected_data_T)
View(countDataMatrix)

# Checking that samples anmes match in both files
all(colnames(countDataMatrix) %in% rownames(selected_meta))
all(colnames(countDataMatrix) == rownames(selected_meta))

# Setting the full design matrix of LRT
design_mat <- model.matrix(~ AGE + TREATMENT + AGE:TREATMENT, selected_meta)

# Convert AGE to factor
selected_meta$AGE <- as.factor(selected_meta$AGE)

# Performing LRT
dds <- DESeqDataSetFromMatrix(countData = countDataMatrix, colData = selected_meta, design = ~ AGE + TREATMENT + AGE:TREATMENT)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ AGE+TREATMENT)

# Getting the LRT result
res_LRT <- results(dds_lrt_time)
res_LRT_df <- as.data.frame(res_LRT)


# Subset the LRT results to return species with padj < 0.05
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="species") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig species that showed a treatment-specific effect at one or more time points after time 0
sigLRT_species <- sig_res_LRT %>% 
  pull(species)

write_xlsx(sig_res_LRT,paste0("../../results/tables/DESeq2_LRT/host_expressed_gene_AGP_sig.xlsx"))






