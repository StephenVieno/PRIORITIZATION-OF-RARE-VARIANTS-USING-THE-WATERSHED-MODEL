
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


# Create gene list of priotized and non-priotized genes 
gene_case <- schema_gene |> 
  janitor::clean_names() |> 
  dplyr::filter(p_meta < 0.10 | p_meta > 0.5) |> 
  mutate(gene_case = ifelse(p_meta < 0.10, "Yes", "No")) 

# Watershed variants within priotized and non-priotized genes  
Schema_GRCh38_watershed <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  inner_join(gene_case, by = c("Gene" = "gene_id"))

# Check difference between VEP and Watershed
Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |>
  dplyr::filter(gene_id != Gene) |> head()
  
# Watershed variants within priotized genes
Schema_GRCh38_watershed_case <- Schema_GRCh38_watershed |> 
  mutate(status = ifelse(gene_case == "Yes", "watershed-case", "No")) |> 
  dplyr::filter(status == "watershed-case") |> 
  distinct(RV, .keep_all = TRUE) |> 
  dplyr::filter(!is.na(est))

# Watershed variants within non-priotized genes
Schema_GRCh38_watershed_control <- Schema_GRCh38_watershed |> 
  mutate(status = ifelse(gene_case == "No", "watershed-control", "No")) |> 
  dplyr::filter(status == "watershed-control") |> 
  distinct(RV, .keep_all = TRUE) |> 
  dplyr::filter(!is.na(est))

# Non-Watershed variants within priotized and non-priotized genes  
Schema_GRCh38_non_watershed <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "No") |> 
  inner_join(gene_case, by = "gene_id") 

# Non-Watershed variants within priotized genes  
Schema_GRCh38_non_watershed_case <- Schema_GRCh38_non_watershed |> 
  mutate(status = ifelse(gene_case == "Yes", "non-watershed-case", "No")) |> 
  dplyr::filter(status == "non-watershed-case") |> 
  distinct(RV, .keep_all = TRUE) |> 
  dplyr::filter(!is.na(est))

# Non-Watershed variants within non-priotized genes  
Schema_GRCh38_non_watershed_control <- Schema_GRCh38_non_watershed |> 
  mutate(status = ifelse(gene_case == "No", "non-watershed-control", "No")) |> 
  dplyr::filter(status == "non-watershed-control") |> 
  distinct(RV, .keep_all = TRUE) |> 
  dplyr::filter(!is.na(est))

# Variant counts 
dim(Schema_GRCh38_watershed_case)
dim(Schema_GRCh38_watershed_control)
dim(Schema_GRCh38_non_watershed_case)
dim(Schema_GRCh38_non_watershed_control)

# Bind datasets 
Schema_case_control_genes <- bind_rows(Schema_GRCh38_non_watershed_case,
                                       Schema_GRCh38_non_watershed_control, 
                                       Schema_GRCh38_watershed_case, 
                                       Schema_GRCh38_watershed_control)

# Wilcoxon signed-rank test
wilcox.test(Schema_GRCh38_watershed_case$rank_est, 
            Schema_GRCh38_watershed_control$rank_est, 
            paired = FALSE)


wilcox.test(Schema_GRCh38_watershed_case$rank_est, 
            Schema_GRCh38_non_watershed_control$rank_est, 
            paired = FALSE)


wilcox.test(Schema_GRCh38_watershed_case$rank_est, 
            Schema_GRCh38_non_watershed_case$rank_est, 
            paired = FALSE)

# Count of variants for each group 
N_Schema_GRCh38_watershed_case <- count(Schema_GRCh38_watershed_case)
N_Schema_GRCh38_watershed_control <- count(Schema_GRCh38_watershed_control)
N_Schema_GRCh38_non_watershed_case <- count(Schema_GRCh38_non_watershed_case)
N_Schema_GRCh38_non_watershed_control <- count(Schema_GRCh38_non_watershed_control)

# Schema_case_control_plot
Schema_case_control_plot <- Schema_case_control_genes |> 
  ggplot(aes(y = rank_est, 
             x = status)) +
  geom_violin(aes(y = rank_est, 
                  x = status)) + 
  geom_boxplot(aes(y = rank_est, 
                   x = status, 
                   fill = status,
                   alpha = 0.5)) + 
  scale_y_continuous( breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  scale_color_manual(values = c("#65CB5EFF","#440154FF","#FDE725FF", "#355E8DFF"), 
                     aesthetics = c("colour", "fill")) + 
  labs(title = "Distribution of Rare Variant Effect Sizes in Schizophrenia-Associated Genes",
       subtitle = "Wilcoxon rank sum test of variant effect size percentiles (p-values reported in figure)", 
       x = "Variant Group",
       y = "Schizophrenia Effect Size Percentile", 
       fill = "Variant Status")  + 
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(face = "bold",  hjust = 1.5, size = 13)) + 
  geom_signif(
    comparisons = list(c("watershed-case", "watershed-control"), 
                       c("watershed-case", "non-watershed-control"),
                       c("watershed-case", "non-watershed-case")),
    textsize = 3, 
    test = "wilcox.test", 
    y_position = c(1, 1.05, 1.1)
  ) + 
  scale_x_discrete(labels = c("watershed-case" = paste0("Watershed case 
gene variants 
(N = ", N_Schema_GRCh38_watershed_case,")"), 
                              "watershed-control" = paste0("Watershed control 
gene variants 
(N = ", N_Schema_GRCh38_watershed_control,")"),
                              "non-watershed-case" = paste0("Non-Watershed case 
gene variants 
(N = ", N_Schema_GRCh38_non_watershed_case,")"),
                              "non-watershed-control" = paste0("Non-Watershed control 
gene variants 
(N = ", N_Schema_GRCh38_non_watershed_control,")")))


# Schema_case_control_plot
Schema_case_control_plot

# Save plot 
ggsave("Plots/Schema_case_control_plot7x7.pdf", 
       plot = Schema_case_control_plot, 
       width = 7, 
       height = 7)

ggsave("Plots/Schema_case_control_plot.pdf", 
       plot = Schema_case_control_plot, 
       width = 7, 
       height = 14)


