
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

# install.packages("remotes")
# remotes::install_github("const-ae/ggsignif")

library(ggsignif)

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
# Watershed vs Non-Watershed Variants -------------------------------------
# -------------------------------------------------------------------------

gene_case <- schema_gene |> 
  janitor::clean_names() |> 
  dplyr::filter(p_meta < 0.10) 


Schema_GRCh38_watershed <- Schema_GRCh38 |> 
  dplyr::filter(gene_id %in% gene_case$gene_id) |> 
  dplyr::filter(!is.na(est))


count1 <- 
  Schema_GRCh38_watershed |> 
  dplyr::filter(Watershed == "Yes") |> 
  count()

count2 <- 
  Schema_GRCh38_watershed |> 
  dplyr::filter(Watershed == "No") |> 
  count()

Schizo_gene_plot <- Schema_GRCh38_watershed |> 
  ggplot(aes(x = Watershed, y = rank_est)) +
  geom_violin() +
  geom_boxplot(aes(fill = Watershed, 
                   alpha = 0.5)) + 
  geom_signif(
    comparisons = list(c("Yes", "No")),
    textsize = 3, 
    test = "wilcox.test") + 
  scale_color_manual(values = c("#355E8DFF", "#FDE725FF"), 
                     aesthetics = c("colour", "fill")) + 
  theme_minimal() + 
    scale_x_discrete(labels = c("Yes" = paste0("Watershed variants 
in genes associated 
with Schizophrenia
(N = ", count1,")"), 
                              "No" = paste0("Non-Watershed variants 
in genes associated 
with Schizophrenia
(N = ", count2,")"))) + 
  labs(title = "Distribution of Rare Variant Effect Sizes 
in Schizophrenia-Associated Genes",
       x = "Watershed Model Prioritization",
       y = "Schizophrenia Effect Size Percentile",
       color = "Dataset") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 12))

# Schizo gene plot  
Schizo_gene_plot

# Save plot 
ggsave("Schizo_gene_plot.pdf", 
       plot = schizo_gene_plot, 
       width = 6, 
       height = 10)










