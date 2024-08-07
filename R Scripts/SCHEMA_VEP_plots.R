
# PACKAGES ----------------------------------------------------------------
library(conflicted)
library(tidyverse)
library(EnsDb.Hsapiens.v79)
library(clusterProfiler)
library(enrichplot)
library(readr)
library(liftOver)
library(viridis)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

# LOAD DATA ---------------------------------------------------------------

schema_data<-
  read_tsv("/Users/stephenvieno/Projects/ROTATION_2/SCHEMA_DATASET/meta_SCHEMA_variant_results.tsv", 
           col_names = FALSE)
schema_gene <- 
  read_tsv("/Users/stephenvieno/Projects/ROTATION_2/SCHEMA_DATASET/SCHEMA_gene_results.tsv", 
           col_names = TRUE)
uv_watershed <- 
  read.csv("/Users/stephenvieno/Projects/ROTATION_2/WATERSHED_MODEL/SupplementaryTable2_PrioritizedRareVariants.csv")

# CHECK DATA --------------------------------------------------------------

# Assign column names of variant data 
col_names <- c("locus", "alleles", "gene_id", "consequence", 
               "hgvsc", "hgvsp", "cadd", "mpc", "polyphen", "group", 
               "ac_case", "ac_ctrl", "an_case", "an_ctrl", "n_denovos",
               "p", "est", "se", "qp", "i2", "in_analysis", "source",  "k")
colnames(schema_data) <- col_names

# Filter for meta analysis data 
schema_data <- schema_data |> dplyr::filter(group == "meta")

# Check number of variants from meta analysis 

table(schema_data$group)
# meta 
# 8859541 

# Check for missingness 
na_counts <- schema_data %>%
  summarize(across(everything(), ~ sum(is.na(.))))

# Print missingness count 
print(na_counts)

# LIFTOVER  ---------------------------------------------------------------

# Convert to GRanges objects
schema_GRanges <- GRanges(schema_data$locus)

# Working directory of chain files 
workdir <- "/Users/stephenvieno/Projects/supporting_files/chain_files/"

## From UCSC Genome Browser
# "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
# /Users/stephenvieno/Projects/supporting_files/chain_files/hg19ToHg38.over.chain
# /Users/stephenvieno/Projects/supporting_files/chain_files/GRCh37_to_GRCh38.chain

## Lift Over ----------------------------------------------------------------

# Import chain files 
liftOverChain <- import.chain(paste0(workdir,"GRCh37_to_GRCh38.chain"))

# Label SNPs 
GenomicRanges::mcols(schema_GRanges)$SNP <- rownames(schema_data)

schema_data$SNP <- rownames(schema_data)

# Liftover files 
b38GRanges <- liftOver(schema_GRanges, liftOverChain)

b38dat <- as.data.frame(unlist(b38GRanges))

# CHECK LIFTOVER FILES ----------------------------------------------------

# Number of variants in liftover file 
dim(b38dat)
# 8833020       6

# Number of variants prior to liftover  
dim(schema_data)
# 8859541      24

table(b38dat$SNP %in% schema_data$SNP)
# TRUE 
# 8833020 

# Number of variants that successfully transferred over during liftover
table(schema_data$SNP %in% b38dat$SNP)
# FALSE    TRUE 
# 26521 8833020 


# MERGE DATASETS AFTER LIFTOVER  ------------------------------------------

# Calculate minor allele frequency (MAF)
# Filter for MAF < 0.01 
schema_data_rare_variants <- schema_data |> 
  mutate(MAF = (ac_case + ac_ctrl) / (an_case + an_ctrl)) |> 
  dplyr::filter(MAF < 0.01)

# Inner join datasets by SNP labels 
# Calculate the percentile of the effect size estimate
# Label variants prioritized by the Watershed model (Yes/No)
Schema_GRCh38 <- b38dat |> 
  inner_join(schema_data_rare_variants, by = "SNP") |> 
  mutate(RV = paste0("chr",seqnames,":",start)) |> 
  left_join(uv_watershed, by = "RV", relationship = "many-to-many") |> 
  mutate(rank_est = percent_rank(est)) |> 
  mutate(Watershed = ifelse(is.na(Gene), "No", "Yes"))


# -------------------------------------------------------------------------
# Variant Effect Predictor Annotations ------------------------------------
# -------------------------------------------------------------------------

# Filter for prioritized Watershed variants
Schema_GRCh38_Watershed <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") 

# VEP vs RNA boxplot  
VEP_plot1 <- 
Schema_GRCh38_Watershed |> 
  ggplot() +
  geom_jitter(aes(x = consequence, y = RNA, alpha = 0.01, color = consequence)) +
  geom_boxplot(aes(x = consequence, y = RNA, alpha = 0.01), outlier.shape = NA) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = viridis(18), aesthetics = c("colour", "fill")) +
  labs(title = "Distribution of RNA Posterior Probability versus Variant Effect Predictor (VEP) Annotations",
       y = "RNA Posterior Probability",
       x = "Variant Effect Predictor (VEP) Consequences",
       color = "VEP Consequences") +
  guides(alpha = "none")

# Plot 1 
VEP_plot1

# Save VEP plot 1 
dev.off()
ggsave(filename = "VEP_plot1.pdf", 
       plot = VEP_plot1, 
       device = "pdf", 
       width = 12, 
       height = 6)


# VEP vs Methylation boxplot  
VEP_plot2 <- 
Schema_GRCh38_Watershed |> 
  ggplot() +
  geom_jitter(aes(x = consequence, y = Methylation, alpha = 0.01, color = consequence)) +
  geom_boxplot(aes(x = consequence, y = Methylation,  alpha = 0.01), outlier.shape = NA) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = viridis(18), aesthetics = c("colour", "fill")) +
  labs(title = "Distribution of Methylation Posterior Probability versus Variant Effect Predictor (VEP) Annotations",
       y = "Methylation Posterior Probability",
       x = "Variant Effect Predictor (VEP) Consequences",
       color = "VEP Consequences") + 
  guides(alpha = "none") + 
  guides(alpha = "none")

# Plot 2 
VEP_plot2

# Save VEP plot 2
dev.off()
ggsave(filename = "VEP_plot2.pdf", 
       plot = VEP_plot2, 
       device = "pdf", 
       width = 12, 
       height = 6)

# VEP vs Splicing boxplot  
VEP_plot3 <- 
Schema_GRCh38_Watershed |> 
  ggplot() +
  geom_jitter(aes(x = consequence, y = Splicing, alpha = 0.01, color = consequence)) +
  geom_boxplot(aes(x = consequence, y = Splicing, alpha = 0.01), outlier.shape = NA) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = viridis(18), aesthetics = c("colour", "fill")) + 
  labs(title = "Distribution of Splicing Posterior Probability versus Variant Effect Predictor (VEP) Annotations",
       y = "Splicing Posterior Probability",
       x = "Variant Effect Predictor (VEP) Consequences",
       color = "VEP Consequences") + 
  guides(alpha = "none")

# Plot 3
VEP_plot3

# Save VEP plot 3
dev.off()
ggsave(filename = "VEP_plot3.pdf", 
       plot = VEP_plot3, 
       device = "pdf", 
       width = 12, 
       height = 6)

# VEP vs Protein boxplot  
VEP_plot4 <- 
Schema_GRCh38_Watershed |> 
  ggplot() +
  geom_jitter(aes(x = consequence, y = Protein, alpha = 0.01, color = consequence)) +
  geom_boxplot(aes(x = consequence, y = Protein,  alpha = 0.01), outlier.shape = NA) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = viridis(18), aesthetics = c("colour", "fill")) + 
  labs(title = "Distribution of Protein Posterior Probability versus Variant Effect Predictor (VEP) Annotations",
       y = "Protein Posterior Probability",
       x = "Variant Effect Predictor (VEP) Consequences",
       color = "VEP Consequences") + 
  guides(alpha = "none")

# Plot 4
VEP_plot4

# Save VEP plot 4
dev.off()
ggsave(filename = "VEP_plot4.pdf", 
       plot = VEP_plot4, 
       device = "pdf", 
       width = 12, 
       height = 6)

