
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
# Logistic Regression with MAF --------------------------------------------
# -------------------------------------------------------------------------

# Watershed variants
Schema_GRCh38_Watershed <- Schema_GRCh38 |> 
  mutate(MAF_scaled = scale(MAF)) |> 
  dplyr::filter(Watershed == "Yes") |> 
  mutate(RNA_binary = ifelse(RNA > 0.5, 1,0)) |> 
  mutate(Methylation_binary = ifelse(Methylation > 0.5, 1,0)) |> 
  mutate(Splicing_binary = ifelse(Splicing > 0.5, 1,0)) |> 
  mutate(Protein_binary = ifelse(Protein > 0.5, 1,0)) 

# RNA Logit Model 
RNA_logit <- glm(RNA_binary ~ cadd + MAF_scaled, 
                 data = Schema_GRCh38_Watershed, 
                 family = "binomial")

RNA_logit_summary_with_MAF <- summary(RNA_logit)

RNA_logit_summary_with_MAF

write.csv(as.data.frame(RNA_logit_summary_with_MAF$coefficients), 
          "RNA_logit_summary.csv")

# Methylation Logit Model 
Methylation_logit <- glm(Methylation_binary ~ cadd + MAF_scaled, 
                         data = Schema_GRCh38_Watershed,
                         family = "binomial")

Methylation_logit_summary_with_MAF <- summary(Methylation_logit)

Methylation_logit_summary_with_MAF

write.csv(as.data.frame(Methylation_logit_summary_with_MAF$coefficients), 
          "Methylation_logit_summary.csv")

# Splicing Logit Model 
Splicing_logit <- glm(Splicing_binary ~ cadd + MAF_scaled, 
                      data = Schema_GRCh38_Watershed,
                      family = "binomial")

Splicing_logit_summary_with_MAF <- summary(Splicing_logit)

Splicing_logit_summary_with_MAF

write.csv(as.data.frame(Splicing_logit_summary_with_MAF$coefficients), 
          "Splicing_logit_summary.csv")

# Protein Logit Model 
Protien_logit <- glm(Protein_binary ~ cadd + MAF_scaled, 
                     data = Schema_GRCh38_Watershed, 
                     family = "binomial")

Protein_logit_summary_with_MAF <- summary(Protien_logit)

Protein_logit_summary_with_MAF

write.csv(as.data.frame(Protein_logit_summary_with_MAF$coefficients), 
          "Protein_logit_summary.csv")

# Hold MAF_scaled at median 

qunatile_MAF <- quantile(Schema_GRCh38_Watershed$MAF_scaled)

qunatile_MAF

predict_data <- 
Schema_GRCh38_Watershed |> 
  mutate(MAF_scaled = median(MAF_scaled)) |> 
  mutate(MAF_scaled = as.matrix(MAF_scaled)) 

# 25% Quantile for MAF
predict_lower <- 
  Schema_GRCh38_Watershed |> 
  mutate(MAF_scaled = -0.12353961) |> 
  mutate(MAF_scaled = as.matrix(MAF_scaled)) 

# 75% Quantile for MAF 
predict_upper <- 
  Schema_GRCh38_Watershed |> 
  mutate(MAF_scaled = 4.65755305) |> 
  mutate(MAF_scaled = as.matrix(MAF_scaled)) 

  

# Logistic regression predictions for plotting 
Schema_GRCh38_Watershed$prediction_RNA <- predict(RNA_logit, 
                                                    newdata = predict_data,  
                                                    type = "response")
Schema_GRCh38_Watershed$prediction_methyl <- predict(Methylation_logit, 
                                                       newdata = predict_data, 
                                                       type = "response")
Schema_GRCh38_Watershed$prediction_splicing <- predict(Splicing_logit, 
                                                         newdata = predict_data, 
                                                         type = "response")
Schema_GRCh38_Watershed$prediction_protein <- predict(Protien_logit, 
                                                        newdata = predict_data, 
                                                        type = "response")


# 25% Quantile Predictions 
Schema_GRCh38_Watershed$preidiction_lower_RNA <- predict(RNA_logit, 
                                                         newdata = predict_lower,  
                                                         type = "response")
Schema_GRCh38_Watershed$preidiction_lower_Methyl <- predict(Methylation_logit, 
                                                            newdata = predict_lower, 
                                                            type = "response")
Schema_GRCh38_Watershed$preidiction_lower_Splicing <- predict(Splicing_logit, 
                                                              newdata = predict_lower, 
                                                              type = "response")
Schema_GRCh38_Watershed$preidiction_lower_Protein <- predict(Protien_logit, 
                                                             newdata = predict_lower, 
                                                             type = "response")



# 75% Quantile Predictions 
Schema_GRCh38_Watershed$preidiction_upper_RNA <- predict(RNA_logit, 
                                                    newdata = predict_upper,  
                                                    type = "response")
Schema_GRCh38_Watershed$preidiction_upper_Methyl <- predict(Methylation_logit, 
                                                       newdata = predict_upper, 
                                                       type = "response")
Schema_GRCh38_Watershed$preidiction_upper_Splicing <- predict(Splicing_logit, 
                                                         newdata = predict_upper, 
                                                         type = "response")
Schema_GRCh38_Watershed$preidiction_upper_Protein <- predict(Protien_logit, 
                                                        newdata = predict_upper, 
                                                        type = "response")


# Logistic plot 
Logistic_plot <- 
  Schema_GRCh38_Watershed |> 
  ggplot() + 
  geom_point(aes(x = cadd, y = RNA, alpha = 0.5)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(aes(y = prediction_RNA, x = cadd), color = "darkred", size = 1) + 
  geom_line(aes(y = preidiction_upper_RNA, x = cadd), color = "red", linetype = "dashed", size = 0.75) + 
  geom_line(aes(y = preidiction_lower_RNA, x = cadd), color = "red", linetype = "dashed", size = 0.75) + 
  labs(title = "RNA Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "RNA Posterior Probability",
       colour = "Prediction") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  Schema_GRCh38_Watershed |> 
  ggplot() + 
  geom_point(aes(x = cadd, y = Methylation, alpha = 0.5)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(aes(y = prediction_methyl, x = cadd), color = "darkred", size = 1) + 
  geom_line(aes(y = preidiction_upper_Methyl, x = cadd), color = "red", linetype = "dashed", size = 0.75) + 
  geom_line(aes(y = preidiction_lower_Methyl, x = cadd), color = "red", linetype = "dashed", size = 0.75) + 
  ylim(0,1) + 
  labs(title = "Methylation Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "Methylation Posterior Probability",
       colour = "Prediction") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  Schema_GRCh38_Watershed|> 
  ggplot() + 
  geom_point(aes(x = cadd, y = Splicing, alpha = 0.5)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(aes(y = prediction_splicing, x = cadd), color = "darkred", size = 1) + 
  geom_line(aes(y = preidiction_upper_Splicing, x = cadd),color = "red", linetype = "dashed", size = 0.75) + 
  geom_line(aes(y = preidiction_lower_Splicing, x = cadd),color = "red", linetype = "dashed", size = 0.75) + 
  labs(title = "Splicing Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "Splicing Posterior Probability",
       colour = "Prediction") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  Schema_GRCh38_Watershed |> 
  ggplot() + 
  geom_point(aes(x = cadd, y = Protein, alpha = 0.5)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_line(aes(y = preidiction_upper_Protein, x = cadd), color = "red", linetype = "dashed", size = 0.75) + 
  geom_line(aes(y = preidiction_lower_Protein, x = cadd), color = "red", linetype = "dashed", size = 0.75) + 
  geom_line(aes(y = prediction_protein, x = cadd), color = "darkred", size = 1) +
  labs(title = "Protein Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "Protein Posterior Probability",
       colour = "Prediction") +
  theme_bw() + 
  theme(legend.position = "none") 


# Logistic_plot
Logistic_plot

# Logistic plot 
Logistic_plot <- 
  Schema_GRCh38_Watershed |> 
  ggplot() + 
  geom_point(aes(x = cadd, y = RNA, alpha = 0.5)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_line(aes(y = prediction_RNA, x = cadd), color = "red", size = 1) + 
  labs(title = "RNA Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "RNA Posterior Probability",
       colour = "Prediction") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  Schema_GRCh38_Watershed |> 
  ggplot() + 
  geom_point(aes(x = cadd, y = Methylation, alpha = 0.5)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_line(aes(y = prediction_methyl, x = cadd), color = "red", size = 1) + 
  ylim(0,1) + 
  labs(title = "Methylation Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "Methylation Posterior Probability",
       colour = "Prediction") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  Schema_GRCh38_Watershed|> 
  ggplot() + 
  geom_point(aes(x = cadd, y = Splicing, alpha = 0.5)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_line(aes(y = prediction_splicing, x = cadd), color = "red", size = 1) + 
  labs(title = "Splicing Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "Splicing Posterior Probability",
       colour = "Prediction") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  Schema_GRCh38_Watershed |> 
  ggplot() + 
  geom_point(aes(x = cadd, y = Protein, alpha = 0.5)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  geom_line(aes(y = prediction_protein, x = cadd), color = "red", size = 1) +
  labs(title = "Protein Signal versus CADD",
       x = "Combined Annotation Dependent Depletion (CADD)",
       y = "Protein Posterior Probability",
       colour = "Prediction") +
  theme_bw() + 
  theme(legend.position = "none") 

# Save pdf 
ggsave(filename = "Logistic_plot.pdf", 
       plot = Logistic_plot, 
       width = 9, 
       height = 6)

# Logistic plot 
Logistic_plot


