
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

# Filter for ultra rare variants 
Schema_GRCh38_Symbols_UV <- Schema_GRCh38_Symbols |> 
  mutate(UV = ac_case + ac_ctrl) |> 
  dplyr::filter(UV <= 10)

# Filter for odds ratio > 1 
Schema_GRCh38_Symbols_UV_filtered <- 
  Schema_GRCh38_Symbols_UV |> 
  mutate(odds_ratio = (ac_case*an_ctrl)/(ac_ctrl*an_case)) |> 
  dplyr::filter(odds_ratio > 1)

# Filter for coding consequences 
uv_data_symbol <- Schema_GRCh38_Symbols_UV_filtered |> 
  dplyr::filter(consequence == "coding_sequence_variant" | 
                  consequence == "missense_variant_mpc_>=3" |
                  consequence == "start_lost" | 
                  consequence == "missense_variant_mpc_2-3" |
                  consequence == "frameshift_variant" | 
                  consequence == "stop_lost" | 
                  consequence == "protein_altering_variant" | 
                  consequence == "splice_acceptor_variant" | 
                  consequence == "splice_region_variant" | 
                  consequence == "splice_donor_variant" |
                  consequence == "stop_gained") 

# Go Enrichment for Watershed gene list 
go_enrich <- enrichGO(gene = unique(uv_data_symbol$gene_id),
                      OrgDb = "org.Hs.eg.db", 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
# Upsetplot 
upsetplot(go_enrich)

# Barplot 
barplot(go_enrich, showCategory=20) + 
  ggtitle("Ultrarare Variants Gene-Set Gene Ontology Enrichment")

# Dotplot 
dotplot <- dotplot(go_enrich, showCategory=20) + 
  ggtitle("Ultrarare Variants Gene-Set Gene Ontology Enrichment") + 
  theme(plot.title = element_text(size = 10)) + 
  xlab("Gene Ratio")

# Dotplot 
dotplot

# Save dotplot
ggsave("Ultrarare_10_dotplot.pdf", 
       plot = dotplot, 
       width = 7, 
       height = 10)
