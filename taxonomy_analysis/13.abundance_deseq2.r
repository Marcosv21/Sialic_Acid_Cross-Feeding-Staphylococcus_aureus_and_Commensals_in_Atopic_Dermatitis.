# =========================================================
# DESeq2 - Differential abundance of MAGs individually
# Comparing groups: AD vs HC
# Using CoverM abundance table
# =========================================================

library(DESeq2)
library(tidyverse)
# Input files
coverm_file <- "/home/marcos/PRJEB59406/table/mag_abundance/mag_abundance_all.tsv"
# Load CoverM table
abund <- read.table(
  coverm_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
# First column = MAG IDs
mag_ids <- abund[,1]
counts <- abund[,-1]
rownames(counts) <- mag_ids
# Clean sample names
clean_names <- colnames(counts)
# Remove CoverM suffix
clean_names <- gsub(".sorted Read Count", "", clean_names, fixed = TRUE)
# Ensure unique names
clean_names <- make.unique(clean_names)
colnames(counts) <- clean_names
counts <- round(as.matrix(counts))
counts <- counts[rowSums(counts) > 0, ]
sample_names <- colnames(counts)
condition <- ifelse(
  grepl("HC_", sample_names),
  "HC",
  "AD"
)
metadata <- data.frame(
  row.names = sample_names,
  condition = condition
)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)
# dds <- DESeq(dds) PRJNA489681
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, sfType = "poscounts")
res <- results(dds,
               contrast = c("condition", "AD", "HC"),
               alpha = 0.05)

# Convert to dataframe
res_df <- as.data.frame(res)

# Add MAG IDs
res_df$MAG <- rownames(res_df)

# Sort by adjusted p-value
res_df <- res_df %>%
  arrange(padj)
# Significant MAGs
# padj < 0.05
sig_res <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05)
# Save outputs
write.table(
  res_df,
  "/home/marcos/PRJEB59406/table/mag_abundance/deseq2_all_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  sig_res,
  "/home/marcos/PRJEB59406/table/mag_abundance/deseq2_significant_MAGs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# Summary
cat("DESeq2 finished successfully\n")
cat("Total MAGs analyzed:", nrow(res_df), "\n")
cat("Significant MAGs:", nrow(sig_res), "\n")