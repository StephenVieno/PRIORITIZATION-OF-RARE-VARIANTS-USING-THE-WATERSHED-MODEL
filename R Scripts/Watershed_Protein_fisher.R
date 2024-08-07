
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

# Create list of all genes in annotations and Watershed model 
all_genes <-  as.data.frame(unique(c(Schema_GRCh38$gene_id, Schema_GRCh38$Gene)))
colnames(all_genes)[1] <- "Symbols"
all_genes <- all_genes |> 
  dplyr::filter(Symbols != "NA")
# -------------------------------------------------------------------------
# Get gene names of annotations 
gene_symbols <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= all_genes$Symbols, 
                                  keytype = "GENEID", 
                                  columns = c("SYMBOL","GENEID"))
# Label gene symbols 
colnames(gene_symbols) <- c("Symbol", "gene_id")

# Combine gene names 
Schema_GRCh38_Symbols <- Schema_GRCh38 |> left_join(gene_symbols, by = "gene_id")

# Watershed symbols 
genelist_watershed <- 
  Schema_GRCh38_Symbols |> 
  dplyr::filter(Watershed == "Yes") |> 
  dplyr::filter(Protein > 0.5) 


#Gene List 1: MIR137 List 
#Gene List 2: pLI09 List 
#Gene List 3: RBFOX2 List 
#Gene List 4: RBFOX3 List 

gene_list1 <- read.table("/Users/stephenvieno/Projects/ROTATION_2/Gene_Lists/mir137.txt")
gene_list2 <- read.table("/Users/stephenvieno/Projects/ROTATION_2/Gene_Lists/pLI09.txt")
gene_list3 <- read.table("/Users/stephenvieno/Projects/ROTATION_2/Gene_Lists/rbfox2.txt")
gene_list4 <- read.table("/Users/stephenvieno/Projects/ROTATION_2/Gene_Lists/rbfox13.txt")

colnames(gene_list1) <- "Symbol"
colnames(gene_list2) <- "Symbol"
colnames(gene_list3) <- "Symbol"
colnames(gene_list4) <- "Symbol"

# -------------------------------------------------------------------------
# GENE LIST 1 -------------------------------------------------------------
# -------------------------------------------------------------------------
genelist_watershed$Symbol
# Total Number of Genes Symbols 
total <- c(Schema_GRCh38_Symbols$Symbol, gene_list1$Symbol) |> unique() |> length()
# Total for Case Status 
total_filtered <- unique(genelist_watershed$Symbol) |> length()

ifelse(unique(gene_list1$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> table()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1 <- ifelse(unique(gene_list1$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> sum()
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2 <- ifelse(unique(gene_list1$Symbol) %in% unique(genelist_watershed$Symbol), 0, 1) |> sum()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2

ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list1$Symbol), 1, 0) |> table()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3 <- ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list1$Symbol), 0, 1) |> sum()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3

unique(genelist_watershed$Symbol) |> length()

# Genes not in produced gene list (MAF < 10) or Shizophrenia genes from Common Variant studies

count4 <- total - (count1 + count2 + count3)

schizogene_and_genelist <- count1
schizogene_and_non_genelist <- count2
non_schizogene_and_genelist <- count3
non_schizogene_non_genelist <- count4

# FISHER TEST -------------------------------------------------------------

matrix_gene <- matrix(c(count1, count3, count2, count4), ncol = 2)
matrix_gene
test1 <- fisher.test(matrix_gene)
test1

# -------------------------------------------------------------------------
# GENE LIST 2 -------------------------------------------------------------
# -------------------------------------------------------------------------

# Total Number of Genes Symbols 
total <- c(Schema_GRCh38_Symbols$Symbol, gene_list2$Symbol) |> unique() |> length()
# Total for Case Status 
total_filtered <- unique(genelist_watershed$Symbol) |> length()

ifelse(unique(gene_list2$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> table()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1 <- ifelse(unique(gene_list2$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> sum()
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2 <- ifelse(unique(gene_list2$Symbol) %in% unique(genelist_watershed$Symbol), 0, 1) |> sum()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2

ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list2$Symbol), 1, 0) |> table()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3 <- ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list2$Symbol), 0, 1) |> sum()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3

unique(genelist_watershed$Symbol) |> length()

# Genes not in produced gene list (MAF < 10) or Shizophrenia genes from Common Variant studies

count4 <- total - (count1 + count2 + count3)

schizogene_and_genelist <- count1
schizogene_and_non_genelist <- count2
non_schizogene_and_genelist <- count3
non_schizogene_non_genelist <- count4

# FISHER TEST -------------------------------------------------------------

matrix_gene <- matrix(c(count1, count3, count2, count4), ncol = 2)
matrix_gene
test2 <- fisher.test(matrix_gene)
test2

# -------------------------------------------------------------------------
# GENE LIST 3 -------------------------------------------------------------
# -------------------------------------------------------------------------

# Total Number of Genes Symbols 
total <- c(Schema_GRCh38_Symbols$Symbol, gene_list3$Symbol) |> unique() |> length()
# Total for Case Status 
total_filtered <- unique(genelist_watershed$Symbol) |> length()

ifelse(unique(gene_list3$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> table()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1 <- ifelse(unique(gene_list3$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> sum()
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2 <- ifelse(unique(gene_list3$Symbol) %in% unique(genelist_watershed$Symbol), 0, 1) |> sum()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2

ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list3$Symbol), 1, 0) |> table()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3 <- ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list3$Symbol), 0, 1) |> sum()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3

unique(genelist_watershed$Symbol) |> length()

# Genes not in produced gene list (MAF < 10) or Shizophrenia genes from Common Variant studies

count4 <- total - (count1 + count2 + count3)

schizogene_and_genelist <- count1
schizogene_and_non_genelist <- count2
non_schizogene_and_genelist <- count3
non_schizogene_non_genelist <- count4

# FISHER TEST -------------------------------------------------------------

matrix_gene <- matrix(c(count1, count3, count2, count4), ncol = 2)
matrix_gene
test3 <- fisher.test(matrix_gene)
test3


# -------------------------------------------------------------------------
# GENE LIST 4 -------------------------------------------------------------
# -------------------------------------------------------------------------

# Total Number of Genes Symbols 
total <- c(Schema_GRCh38_Symbols$Symbol, gene_list4$Symbol) |> unique() |> length()
# Total for Case Status 
total_filtered <- unique(genelist_watershed$Symbol) |> length()

ifelse(unique(gene_list4$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> table()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1 <- ifelse(unique(gene_list4$Symbol) %in% unique(genelist_watershed$Symbol), 1, 0) |> sum()
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2 <- ifelse(unique(gene_list4$Symbol) %in% unique(genelist_watershed$Symbol), 0, 1) |> sum()
# Number of Shizophrenia genes from Common Variant studies in produced gene list (MAF < 10)
count1
# Number of Shizophrenia genes from Common Variant studies not in produced gene list (MAF < 10)
count2

ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list4$Symbol), 1, 0) |> table()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3 <- ifelse(unique(genelist_watershed$Symbol) %in% unique(gene_list4$Symbol), 0, 1) |> sum()
# Number of genes from produced gene list (MAF < 10) not in Shizophrenia genes from Common Variant studies
count3

unique(genelist_watershed$Symbol) |> length()

# Genes not in produced gene list (MAF < 10) or Shizophrenia genes from Common Variant studies

count4 <- total - (count1 + count2 + count3)

schizogene_and_genelist <- count1
schizogene_and_non_genelist <- count2
non_schizogene_and_genelist <- count3
non_schizogene_non_genelist <- count4

# FISHER TEST -------------------------------------------------------------

matrix_gene <- matrix(c(count1, count3, count2, count4), ncol = 2)
matrix_gene
test4 <- fisher.test(matrix_gene)
test4

# -------------------------------------------------------------------------

# Create summary data frame 
summary_fisher <- data.frame(
  gene_list = paste0("gene list ", seq(1:4)), 
  list = c("MIR137 List", "pLI09 List", "RBFOX2 List", "RBFOX3 List"),
  odds_ratio = c(test1$estimate[1], test2$estimate[1], test3$estimate[1], test4$estimate[1]), 
  p_value = c(test1$p.value[1], test2$p.value[1], test3$p.value[1], test4$p.value[1]), 
  lower_CI = c(test1$conf.int[1], test2$conf.int[1], test3$conf.int[1], test4$conf.int[1]),
  upper_CI = c(test1$conf.int[2], test2$conf.int[2], test3$conf.int[2], test4$conf.int[2])
)

head(summary_fisher)

# Multiple testing correction 
summary_fisher$BH_p_values <- p.adjust(summary_fisher$p_value, method = "BH")

# Save data frame
write.csv(summary_fisher, 
          "Watershed_Protein_fisher.csv", 
          row.names = FALSE)

