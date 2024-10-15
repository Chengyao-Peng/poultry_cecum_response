### Bioconductor and CRAN libraries used
library(ggfortify)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(readr)
library(readxl)
library(microbiome)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(gtools)
library(writexl)

# Define age, measure type and cutoff
padj.cutoff <- 0.05
lfc.cutoff <- 0
day <- 3
measure_type <- "Mucosa"
gut_section <- "Cecum"


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

num_samples <- nrow(na.omit(selected_meta))


# Read count table
data <- read_tsv("../../data/processed/rna_seq_filtered_low.tsv")
data <- data %>% column_to_rownames(var="sample_id")
View(data)

# Match the count data and the metadata
selected_data <- data[rownames(data) %in% rownames(selected_meta), ]
selected_data <- selected_data[match(rownames(selected_meta), rownames(selected_data)), ]

# Convert the format of the count data
selected_data_T <- t(selected_data)
countDataMatrix <- as.matrix(selected_data_T)

# checking that samples anmes match in both files
all(colnames(countDataMatrix) %in% rownames(selected_meta))
all(colnames(countDataMatrix) == rownames(selected_meta))

# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = countDataMatrix, colData = selected_meta, design = ~ TREATMENT)
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

#  **Optional step** - Output normalized counts to save as a file to access outside RStudio
normalized_counts <- counts(dds, normalized=TRUE)

# Plot dispersion estimates
plotDispEsts(dds)

# Output results of Wald test for contrast
contrast_AGP_CTR <- c("TREATMENT", "AGP", "CTR")
res_AGP_CTR <- results(dds, contrast = contrast_AGP_CTR)
res_AGP_CTR <- lfcShrink(dds, coef = "TREATMENT_AGP_vs_CTR", res=res_AGP_CTR, type="apeglm")
summary(res_AGP_CTR)

# Turn the results object into a data frame and output
res_AGP_CTR_df <- data.frame(res_AGP_CTR)
AGP_CTR_output_all <- paste0("../../results/tables/RNA_seq/", gut_section, "_", measure_type, "_", day, "_AGP_CTR_all.xlsx")
res_AGP_CTR_df1 <- cbind(" "=rownames(res_AGP_CTR_df), res_AGP_CTR_df)
write_xlsx(res_AGP_CTR_df1, AGP_CTR_output_all)

# Subset the significant results
sig_res_AGP_CTR <- filter(res_AGP_CTR_df, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
AGP_CTR_output_sig <- paste0("../../results/tables/RNA_seq", gut_section, "_", measure_type, "_", day, "_AGP_CTR_shrunk_stat_results.xlsx")
sig_res_AGP_CTR1 <- cbind(" "=rownames(sig_res_AGP_CTR), sig_res_AGP_CTR)
write_xlsx(sig_res_AGP_CTR1,AGP_CTR_output_sig)
sig_num_AGP_CTR = nrow(sig_res_AGP_CTR)
sig_num_AGP_CTR

# -------------------------------------------------------------------------
