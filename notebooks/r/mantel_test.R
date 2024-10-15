# Load necessary libraries using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, RColorBrewer, vegan, pheatmap)

# Helper function to read and format data
read_and_format <- function(file_path) {
  df <- read.delim(file_path)
  rownames(df) <- df[, 1]
  df <- df[, -1]
  return(df)
}

# File paths
mtg_digesta_path <- "../../data/processed/metagenomics_digesta_species_samples.tsv"
mtg_mucosa_path <- "../../data/processed/metagenomics_mucosa_species_samples.tsv"
host_rna_path <- "../../data/processed/host_RNA_samples.tsv"

# Load the datasets
mtg_digesta.sp <- read_and_format(mtg_digesta_path)
mtg_mucosa.sp <- read_and_format(mtg_mucosa_path)
host.rna <- read_and_format(host_rna_path)

# Set names and count features
name <- c("Digesta microbiota species", "Mucosa microbiota species", "Host transcripts")
num_features <- c(nrow(mtg_digesta.sp), nrow(mtg_mucosa.sp), nrow(host.rna))

# Put all datasets into a list
dat <- list(mtg_digesta.sp, mtg_mucosa.sp, host.rna)

# Check for NA values
map(dat, ~sum(is.na(.)))

# Sort rows and columns by column names
dat <- purrr::map(dat, ~ .[, order(colnames(.))])

# Calculate Spearman correlation distance
cor.spe <- purrr::map(dat, ~ 1 - cor(.x, method = "spearman"))

# Create pairs for all combinations
pair <- combn(1:length(dat), 2)

# Initialize an empty list for Mantel test results
mantel.spe <- vector("list", ncol(pair))

# Perform Mantel tests for all pairs
for (j in 1:ncol(pair)) {
  i <- pair[1, j]
  ii <- pair[2, j]
  mantel.spe[[j]] <- mantel(cor.spe[[i]], cor.spe[[ii]])$statistic
}

# Initialize matrices for Mantel distances
mat.spe <- matrix(1, nrow = length(dat), ncol = length(dat))
colnames(mat.spe) <- name
rownames(mat.spe) <- name

# Fill in Mantel distances
for (j in 1:ncol(pair)) {
  i <- pair[1, j]
  ii <- pair[2, j]
  mat.spe[i, ii] <- mantel.spe[[j]]
  mat.spe[ii, i] <- mantel.spe[[j]]
}

view(mat.spe)
