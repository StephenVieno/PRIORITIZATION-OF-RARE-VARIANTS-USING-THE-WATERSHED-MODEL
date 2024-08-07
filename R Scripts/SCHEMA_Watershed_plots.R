
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
# Watershed Model Variants  -----------------------------------------------
# -------------------------------------------------------------------------

# Prioritized variants with Posterior Probability > 0.5 
Schema_GRCh38_Watershed_0.5 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  dplyr::filter(RNA > 0.5 | Methylation > 0.5 | Splicing > 0.5| Protein > 0.5) |> 
  distinct(RV, .keep_all = TRUE) |> 
  mutate(Prioritized = ifelse(Watershed == "Yes", "Watershed_0.5", "No")) 

# Prioritized variants with Posterior Probability > 0.9 
Schema_GRCh38_Watershed_0.9 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  dplyr::filter(RNA > 0.9 | Methylation > 0.9 | Splicing > 0.9| Protein > 0.9) |> 
  distinct(RV, .keep_all = TRUE) |> 
  mutate(Prioritized = ifelse(Watershed == "Yes", "Watershed_0.9", "No")) 

# Prioritized variants with RNA > 0.5 
Schema_GRCh38_RNA_0.5 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(RNA > 0.5, "RNA_0.5", "No")) |> 
  dplyr::filter(Prioritized == "RNA_0.5") |> 
  distinct(RV, .keep_all = TRUE) 

# Prioritized variants with Methylation > 0.5 
Schema_GRCh38_Methyl_0.5 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Methylation > 0.5, "Methylation_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Methylation_0.5") |> 
  distinct(RV, .keep_all = TRUE) 

# Prioritized variants with Splicing > 0.5 
Schema_GRCh38_Splicing_0.5 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Splicing > 0.5, "Splicing_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Splicing_0.5") |> 
  distinct(RV, .keep_all = TRUE) 

# Prioritized variants with Protein > 0.5 
Schema_GRCh38_Protein_0.5 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Protein > 0.5, "Protein_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Protein_0.5") |> 
  distinct(RV, .keep_all = TRUE)

# Prioritized variants with RNA > 0.9 
Schema_GRCh38_RNA_0.9 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(RNA > 0.9, "RNA_0.9", "No")) |> 
  dplyr::filter(Prioritized == "RNA_0.9") |> 
  distinct(RV, .keep_all = TRUE) 

# Prioritized variants with Methylation > 0.9 
Schema_GRCh38_Methyl_0.9 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Methylation > 0.9, "Methylation_0.9", "No")) |> 
  dplyr::filter(Prioritized == "Methylation_0.9") |> 
  distinct(RV, .keep_all = TRUE) 

# Prioritized variants with Splicing > 0.9 
Schema_GRCh38_Splicing_0.9 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Splicing > 0.9, "Splicing_0.9", "No")) |> 
  dplyr::filter(Prioritized == "Splicing_0.9") |> 
  distinct(RV, .keep_all = TRUE) 

# Prioritized variants with Protein > 0.9 
Schema_GRCh38_Protein_0.9 <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Protein > 0.9, "Protein_0.9", "No")) |> 
  dplyr::filter(Prioritized == "Protein_0.9") |> 
  distinct(RV, .keep_all = TRUE) 

# Non-prioritized variants
Schema_GRCh38_non_watershed <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "No") |> 
  dplyr::filter(!is.na(est))

dim(Schema_GRCh38_non_watershed)
head(Schema_GRCh38_non_watershed)

# Check that rank_est matches est 
Schema_GRCh38 |> 
  dplyr::filter(Watershed == "No") |> 
  dplyr::filter(!is.na(rank_est)) |> 
  dim()

# Bind prioritized variants together 
Prioritized_data <- bind_rows(Schema_GRCh38_Watershed_0.5, 
                              Schema_GRCh38_Watershed_0.9, 
                              Schema_GRCh38_RNA_0.5, 
                              Schema_GRCh38_Methyl_0.5, 
                              Schema_GRCh38_Splicing_0.5, 
                              Schema_GRCh38_Protein_0.5,
                              Schema_GRCh38_RNA_0.9, 
                              Schema_GRCh38_Methyl_0.9, 
                              Schema_GRCh38_Splicing_0.9, 
                              Schema_GRCh38_Protein_0.9)

# Number unique of prioritized variants
number0.5 <- unique(Schema_GRCh38_Watershed_0.5$RV) |> length()
number0.9 <- unique(Schema_GRCh38_Watershed_0.9$RV) |> length()
number0.5_RNA <- unique(Schema_GRCh38_RNA_0.5$RV) |> length()
number0.5_Methyl <- unique(Schema_GRCh38_Methyl_0.5$RV) |> length()
number0.5_Splicing <- unique(Schema_GRCh38_Splicing_0.5$RV) |> length()
number0.5_Protein <- unique(Schema_GRCh38_Protein_0.5$RV) |> length()
number0.9_RNA <- unique(Schema_GRCh38_RNA_0.9$RV) |> length()
number0.9_Methyl <- unique(Schema_GRCh38_Methyl_0.9$RV) |> length()
number0.9_Splicing <- unique(Schema_GRCh38_Splicing_0.9$RV) |> length()
number0.9_Protein <- unique(Schema_GRCh38_Protein_0.9$RV) |> length()

# Wilcoxon Rank Sum
wilcox.test(Schema_GRCh38_Watershed_0.5$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_Watershed_0.9$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_RNA_0.5$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_Splicing_0.5$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_Protein_0.5$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_Watershed_0.5$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_RNA_0.9$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_Methyl_0.9$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_Splicing_0.9$est, 
            Schema_GRCh38_non_watershed$est, 
            paired = FALSE)

wilcox.test(Schema_GRCh38_Protein_0.9$rank_est, 
            Schema_GRCh38_non_watershed$rank_est, 
            paired = FALSE)

# Watershed plot 
watershed_plot <- Prioritized_data |> 
  ggplot() + 
  geom_violin(aes(y = rank_est, 
                  x = Prioritized, 
                  alpha = 0.5)) + 
  geom_boxplot(aes(y = rank_est, 
                   x = Prioritized, 
                   fill = Prioritized, 
                   alpha = 0.5)) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") + 
  scale_color_manual(values = c("#440154FF", "#481769FF", 
                                "#3D4E8AFF", "#355E8DFF", 
                                "#23898EFF", "#1F978BFF",
                                "#46C06FFF", "#65CB5EFF",
                                "#D8E219FF", "#FDE725FF"), 
                     aesthetics = c("colour", "fill")) + 
  labs(title = "Distribution of Rare Variant Effect Sizes Prioritized by the Watershed Model at Different Posterior Probabilities",
       x = "Posterior Probabilities",
       y = "Schizophrenia Effect Size Percentile",
       color = "Dataset") +
  guides(alpha = FALSE, 
         fill = FALSE) +
  scale_x_discrete(labels = c("Watershed_0.5" = paste0("Watershed > 0.5 (N = ", number0.5,")"), 
                              "Watershed_0.9" = paste0("Watershed > 0.9 (N = ", number0.9,")"),
                              "RNA_0.5" = paste0("RNA > 0.5 (N = ", number0.5_RNA,")"),
                              "Methylation_0.5" = paste0("Methylation > 0.5 (N = ", number0.5_Methyl,")"),
                              "Splicing_0.5" = paste0("Splicing > 0.5 (N = ", number0.5_Splicing,")"),
                              "Protein_0.5" = paste0("Protein > 0.5 (N = ", number0.5_Protein,")"),
                              "RNA_0.9" = paste0("RNA > 0.9 (N = ", number0.9_RNA,")"),
                              "Methylation_0.9" = paste0("Methylation > 0.9 (N = ", number0.9_Methyl,")"),
                              "Splicing_0.9" = paste0("Splicing > 0.9 (N = ", number0.9_Splicing,")"),
                              "Protein_0.9" = paste0("Protein > 0.9 (N = ", number0.9_Protein,")")))

# Watershed plot 
watershed_plot

# Save plot 
ggsave("Watershed_plot.pdf", 
       plot = watershed_plot, 
       width = 10, 
       height = 7)


