
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
# Watershed in Case and Control Genes -------------------------------------
# -------------------------------------------------------------------------

# Generate gene lists 
gene_case <- schema_gene |> 
  janitor::clean_names() |> 
  dplyr::filter(p_meta < 0.10 | p_meta > 0.5) |> 
  mutate(gene_case = ifelse(p_meta < 0.10, "Yes", "No")) 

# Watershed variants within priotized and non-priotized genes  
Schema_GRCh38_watershed <- Schema_GRCh38 |> 
  dplyr::filter(Watershed == "Yes") |> 
  inner_join(gene_case, by = c("Gene" = "gene_id"))

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

# -------------------------------------------------------------------------

# Variants prioritized by the watershed model at PP > 0.5 in case genes
Schema_GRCh38_Watershed_0.5_case <- Schema_GRCh38_watershed_case |> 
  dplyr::filter(RNA > 0.5 | Methylation > 0.5 | Splicing > 0.5| Protein > 0.5) |> 
  mutate(Prioritized = ifelse(Watershed == "Yes", "Watershed_0.5", "No")) |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at RNA PP > 0.5 in case genes
Schema_GRCh38_RNA_0.5_case <- Schema_GRCh38_watershed_case |> 
  mutate(Prioritized = ifelse(RNA > 0.5, "RNA_0.5", "No")) |> 
  dplyr::filter(Prioritized == "RNA_0.5") |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at Methylation PP > 0.5 in case genes
Schema_GRCh38_Methyl_0.5_case <- Schema_GRCh38_watershed_case |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Methylation > 0.5, "Methylation_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Methylation_0.5") |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at PP > Splicing 0.5 in case genes
Schema_GRCh38_Splicing_0.5_case <- Schema_GRCh38_watershed_case |> 
  mutate(Prioritized = ifelse(Splicing > 0.5, "Splicing_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Splicing_0.5") |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at Protein PP > 0.5 in case genes
Schema_GRCh38_Protein_0.5_case <- Schema_GRCh38_watershed_case |> 
  mutate(Prioritized = ifelse(Protein > 0.5, "Protein_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Protein_0.5") |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at PP > 0.5 in control genes
Schema_GRCh38_Watershed_0.5_control <- Schema_GRCh38_watershed_control |> 
  dplyr::filter(RNA > 0.5 | Methylation > 0.5 | Splicing > 0.5| Protein > 0.5) |> 
  mutate(Prioritized = ifelse(Watershed == "Yes", "Watershed_0.5", "No")) |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at RNA PP > 0.5 in control genes
Schema_GRCh38_RNA_0.5_control <- Schema_GRCh38_watershed_control |> 
  mutate(Prioritized = ifelse(RNA > 0.5, "RNA_0.5", "No")) |> 
  dplyr::filter(Prioritized == "RNA_0.5") |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at Methylation PP > 0.5 in control genes
Schema_GRCh38_Methyl_0.5_control <- Schema_GRCh38_watershed_control |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(Prioritized = ifelse(Methylation > 0.5, "Methylation_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Methylation_0.5") |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at Splicing PP > 0.5 in control genes
Schema_GRCh38_Splicing_0.5_control <- Schema_GRCh38_watershed_control |> 
  mutate(Prioritized = ifelse(Splicing > 0.5, "Splicing_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Splicing_0.5") |> 
  dplyr::filter(est != "NA") 

# Variants prioritized by the watershed model at Protein PP > 0.5 in control genes
Schema_GRCh38_Protein_0.5_control <- Schema_GRCh38_watershed_control |> 
  mutate(Prioritized = ifelse(Protein > 0.5, "Protein_0.5", "No")) |> 
  dplyr::filter(Prioritized == "Protein_0.5") |> 
  dplyr::filter(est != "NA") 

# -------------------------------------------------------------------------

# Bind together variant data sets for plotting 
Prioritized_data_0.5 <- bind_rows(Schema_GRCh38_Watershed_0.5_case, 
                                  Schema_GRCh38_Watershed_0.5_control,
                                  Schema_GRCh38_RNA_0.5_case, 
                                  Schema_GRCh38_Methyl_0.5_case, 
                                  Schema_GRCh38_Splicing_0.5_case, 
                                  Schema_GRCh38_Protein_0.5_case,
                                  Schema_GRCh38_RNA_0.5_control, 
                                  Schema_GRCh38_Methyl_0.5_control, 
                                  Schema_GRCh38_Splicing_0.5_control, 
                                  Schema_GRCh38_Protein_0.5_control)

# Create order for plotting 
order <- c("watershed-control&Methylation_0.5", "watershed-case&Methylation_0.5", 
           "watershed-control&Protein_0.5", "watershed-case&Protein_0.5",
           "watershed-control&RNA_0.5", "watershed-case&RNA_0.5", 
           "watershed-control&Splicing_0.5", "watershed-case&Splicing_0.5",
           "watershed-control&Watershed_0.5", "watershed-case&Watershed_0.5")

# Comparisons for statistic tests 
my_comparisons <- 
  list(c("watershed-control&Methylation_0.5", "watershed-case&Methylation_0.5"), 
       c("watershed-control&Protein_0.5", "watershed-case&Protein_0.5"),
       c("watershed-control&RNA_0.5", "watershed-case&RNA_0.5"), 
       c("watershed-control&Splicing_0.5", "watershed-case&Splicing_0.5"), 
       c("watershed-control&Watershed_0.5", "watershed-case&Watershed_0.5"))

# Create data frame for plotting 
plot_data <- Prioritized_data_0.5 |> 
  mutate(status_Prioritized = paste0(status,"&",Prioritized)) |> 
  mutate(status_Prioritized = factor(status_Prioritized, levels = order))  

# Count of variants for each status definition 
count1 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-control&Methylation_0.5") |> 
  count()

count2 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-case&Methylation_0.5") |> 
  count()

count3 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-control&Protein_0.5") |> 
  count()

count4 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-case&Protein_0.5") |> 
  count()

count5 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-control&RNA_0.5") |> 
  count()

count6 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-case&RNA_0.5") |> 
  count()

count7 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-control&Splicing_0.5") |> 
  count()

count8 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-case&Splicing_0.5") |> 
  count()

count9 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-control&Watershed_0.5") |> 
  count()

count10 <- plot_data |> 
  dplyr::filter(status_Prioritized == "watershed-case&Watershed_0.5") |> 
  count()

# Schema_variants_0.5_plot
Schema_variants_0.5_plot <- plot_data |> 
  mutate(status_Prioritized = paste0(status,"&",Prioritized)) |> 
  mutate(status_Prioritized = factor(status_Prioritized, levels = order)) |> 
  ggplot(aes(y = rank_est, x = status_Prioritized)) + 
  geom_boxplot(aes( 
                   fill = status, 
                   alpha = 0.5)) + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 10)) + 
  scale_color_manual(values = c("#FDE725FF", "#355E8DFF"), 
                     aesthetics = c("colour", "fill")) + 
  labs(title = "Distribution of Rare Variant Effect Sizes Prioritized by the Watershed Model in Schizophrenia-Associated Genes at Posterior Probabilities > 0.5",
       subtitle = "Wilcoxon rank sum test of variant effect size percentiles (p-values reported in figure)", 
       x = "",
       y = "Schizophrenia Effect Size Percentile",
       color = "Dataset") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(face = "bold")) + 
  geom_signif(
    comparisons = my_comparisons,
    textsize = 3, 
    test = "wilcox.test") + 
  scale_x_discrete(labels = c("watershed-control&Methylation_0.5" = paste0("Control genes
& Methylation > 0.5
(N = ", count1,")"), 
                              "watershed-case&Methylation_0.5" = paste0("Case genes
& Methylation > 0.5
(N = ", count2,")"),
                              "watershed-control&Protein_0.5" = paste0("Control genes
& Protein > 0.5
(N = ", count3,")"),
                              "watershed-case&Protein_0.5" = paste0("Case genes
& Protein > 0.5
(N = ", count4,")"),
                              "watershed-control&RNA_0.5" = paste0("Control genes
& RNA > 0.5
(N = ", count5,")"),
                              "watershed-case&RNA_0.5" = paste0("Case genes
& RNA > 0.5
(N = ", count6,")"),
                              "watershed-control&Splicing_0.5" = paste0("Control genes
& Splicing > 0.5
(N = ", count7,")"),
                              "watershed-case&Splicing_0.5" = paste0("Case genes
& Splicing > 0.5
(N = ", count8,")"),
                              "watershed-control&Watershed_0.5" = paste0("Control genes
& Watershed > 0.5
(N = ", count9,")"),
                              "watershed-case&Watershed_0.5" = paste0("Case genes
& Watershed > 0.5
(N = ", count10,")")))
                   
                              
#  Schema_variants_0.5_plot
Schema_variants_0.5_plot

# Save plot 
ggsave("Plots/Schema_variants_0.5_plot.pdf", 
       plot = Schema_variants_0.5_plot, 
       width = 10, 
       height = 7)



# Multiple testing correction 
test1 <- 
  wilcox.test(Schema_GRCh38_Watershed_0.5_case$rank_est, 
              Schema_GRCh38_Watershed_0.5_control$rank_est)

test2 <- 
  wilcox.test(Schema_GRCh38_RNA_0.5_case$rank_est, 
              Schema_GRCh38_RNA_0.5_control$rank_est)

test3 <- 
  wilcox.test(Schema_GRCh38_Methyl_0.5_case$rank_est, 
              Schema_GRCh38_Methyl_0.5_control$rank_est)

test4 <- 
  wilcox.test(Schema_GRCh38_Splicing_0.5_case$rank_est, 
              Schema_GRCh38_Splicing_0.5_control$rank_est)

test5 <- 
  wilcox.test(Schema_GRCh38_Protein_0.5_case$rank_est, 
              Schema_GRCh38_Protein_0.5_control$rank_est)

tests <- c(test1$p.value, 
  test2$p.value,
  test3$p.value,
  test4$p.value,
  test5$p.value)

adjusted_p_values <- p.adjust(tests, method = "BH")

print(adjusted_p_values)
#  0.48396251 0.06247799 0.97215367 0.48396251 0.29310209

# -------------------------------------------------------------------------





