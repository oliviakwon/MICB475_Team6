# ##### Install packages #####
# # Start by installing all necessary packages when asked if you want to install
# # from source, please just type Yes in the terminal below
# 
# # If you don't have BiocManager, here is the code to install it
# # A lot of you probably already have this so you can skip
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # Create a list of all the packages you need to install
# pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
#           "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
#           "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")
# 
# # Use the above list to install all the packages using a for loop
# for (pkg in pkgs) {
#   if (!requireNamespace(pkg, quietly = TRUE))
#     BiocManager::install(pkg)
# }
# # when asked if you want to update all, some or none please type "n" for none
# 
# # After installing all of its above dependencies, install ggpicrust2
# install.packages("ggpicrust2")
# install.packages("readxl")

#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(ggrepel)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(readxl)
library(patchwork)
library(ggplot2)
library(tidyr)
library(scales)
library(phyloseq)
library(ggsignif)
library(ggh4x)
library(ggside)
library(reshape2)

# Load phyloseq object from Aim 2
load("FD_final_only_separate.Rdata")

#duplicate phylosesq object
phyloseq_FD<- FD_final_only_separate 

#extract the sample data
sample_phyloseq_FD <- sample_data(phyloseq_FD)
sample_phyloseq_FD <- as.data.frame(sample_phyloseq_FD)
sample_phyloseq_FD$'sample-id' <- rownames(sample_phyloseq_FD)

#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "picrust2_out_pipeline/pathabun_exported/pathway_abundance.tsv"
#Skipped first row (passage: exported from biome file)
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
abundance_data = as.data.frame(abundance_data)

#Import metadata file, no need to filter yet
metadata <- read_excel("dysautonomia_metadata.xlsx")
metadata <- metadata[metadata$`sample-id` %in% sample_phyloseq_FD$`sample-id`, ]

#Looking at Sex

#Remove NAs for your column of interest in this case subject
metadata<-metadata[!is.na(metadata$Sex),]

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = metadata$'sample-id'
sample_names = append(sample_names, "#OTU ID")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata = metadata[metadata$`sample-id` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
# abundance_data_filtered[ , 'pathway'] = NA

abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
                                        metadata = metadata, group = "Sex", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Make sure all p_adjust_values are in numeric form (some are in String form, for ex. "9.40629626379781e-05")
abundance_daa_results_df$p_adjust <- as.numeric(abundance_daa_results_df$p_adjust)

# Filter p-values to only significant ones
feature_with_p_adjust_0.05 <- abundance_daa_results_df %>% filter(p_adjust < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_adjust_0.05, metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_adjust_0.05)

colnames(abundance_data_filtered)[colnames(abundance_data_filtered) == "#OTU ID"] <- "pathway"

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_adjust_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
num_columns <- ncol(abundance_desc)
print(num_columns) #67 so 67 - 6("feature, "method", "group1", "group2", "p_values", "adj_method", "p_adjust", "description") = 61
abundance_desc = abundance_desc[,-c(61:ncol(abundance_desc))] 

# Generate a heatmap
"feature" %in% colnames(abundance_desc) #check
nrow(abundance_desc) > 0 #check

# Combine 'Sex' and 'Genotype' into a new column 'Group'
sex_gen_metadata <- metadata
sex_gen_metadata$Group <- paste(metadata$Sex, metadata$Genotype, sep = "_")

pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), 
                metadata = sex_gen_metadata, 
                group = "Group")

### Generate pathway PCA plot for mutants

control_metadata <- metadata[metadata$Genotype == 'Control',]
mutant_metadata <- metadata[metadata$Genotype == 'Mutant',]

#Filtering the abundance table to only include samples that are in the filtered metadata
mutant_sample_names = mutant_metadata$'sample-id'
mutant_sample_names = append(mutant_sample_names, "#OTU ID")
mutant_abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% mutant_sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
mutant_abundance_data_filtered =  mutant_abundance_data_filtered[, colSums(mutant_abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(mutant_abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
mutant_abun_samples = rownames(t(mutant_abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
mutant_metadata = mutant_metadata[mutant_metadata$`sample-id` %in% mutant_abun_samples,] #making sure the filtered metadata only includes these samples

colnames(mutant_abundance_data_filtered)[colnames(mutant_abundance_data_filtered) == "#OTU ID"] <- "pathway"
mutant_abundance_data_filtered <- mutant_abundance_data_filtered %>% column_to_rownames("pathway")

# fo r re-scale error
# Calculate row standard deviations
row_sds <- apply(mutant_abundance_data_filtered, 1, sd)

# Filter out rows with zero standard deviation
filtered_data <- mutant_abundance_data_filtered[row_sds != 0, ]
pathway_pca(abundance = filtered_data, metadata = mutant_metadata, group = "Sex")

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
mutant_abundance_data_filtered$pathway = rownames(mutant_abundance_data_filtered)
rownames(mutant_abundance_data_filtered) = NULL
mutant_abundance_daa_results_df <- pathway_daa(abundance = mutant_abundance_data_filtered %>% column_to_rownames("pathway"),
                                               metadata = mutant_metadata, group = "Sex", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
mutant_metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                              daa_results_df = mutant_abundance_daa_results_df, ko_to_kegg = FALSE)

mutant_abundance_daa_results_df$p_adjust <- as.numeric(mutant_abundance_daa_results_df$p_adjust)
# Filter p-values to only significant ones
mutant_feature_with_p_adjust_0.05 <- mutant_abundance_daa_results_df %>% filter(p_adjust < 0.05)

#Changing the pathway column to description for the results 
mutant_feature_desc = inner_join(mutant_feature_with_p_adjust_0.05,mutant_metacyc_daa_annotated_results_df, by = "feature")
mutant_feature_desc$feature = mutant_feature_desc$description
mutant_feature_desc = mutant_feature_desc[,c(1:7)]
colnames(mutant_feature_desc) = colnames(mutant_feature_with_p_adjust_0.05)

### Generating a bar plot representing log2FC from the custom deseq2 function
# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")

# Running DEseq2 function
res =  DEseq2_function(mutant_abundance_data_filtered, mutant_metadata, "Sex")
res$feature =rownames(res)
res_desc = inner_join(res,mutant_metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# filter to only include significant pathways by Log2fold change
# keep only greater than 2 fold change or less than -2
res_desc$padj <- as.numeric(res_desc$padj)
res_desc$log2FoldChange <- as.numeric(res_desc$log2FoldChange)
sig_res = res_desc %>%
  filter(padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2 ))

### Mutant Heatmap with log2FC cutoff & p < 0.05
#Changing the pathway column to description for the abundance table
mutant_abundance = mutant_abundance_data_filtered %>% filter(pathway %in% mutant_feature_with_p_adjust_0.05$feature)
mutant_abundance <- mutant_abundance[, c("pathway", setdiff(names(mutant_abundance), "pathway"))]
colnames(mutant_abundance)[colnames(mutant_abundance) == "pathway"] <- "feature"
mutant_abundance_desc = inner_join(mutant_abundance,mutant_metacyc_daa_annotated_results_df, by = "feature")
mutant_abundance_desc$feature = mutant_abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
num_columns <- ncol(mutant_abundance_desc)
print(num_columns) #42 so 42 - 6("method", "group1", "group2", "p_values", "adj_method", "p_adjust", "description") = 36
mutant_abundance_desc = mutant_abundance_desc[,-c(36:ncol(mutant_abundance_desc))] 

colnames(sig_res)[colnames(sig_res) == "feature"] <- "old_feature"
colnames(sig_res)[colnames(sig_res) == "description"] <- "feature"

log_cutoff_abundance_desc <- mutant_abundance_desc[mutant_abundance_desc$feature %in% sig_res$feature, ]

rownames(log_cutoff_abundance_desc) = NULL
# Generate a heatmap for mutant genotype & Sex
"feature" %in% colnames(log_cutoff_abundance_desc) #check
nrow(log_cutoff_abundance_desc) > 0 #check
log_cutoff_abundance_desc <- log_cutoff_abundance_desc %>% column_to_rownames("feature")

pathway_heatmap(abundance = log_cutoff_abundance_desc, 
                metadata = mutant_metadata, 
                group = "Sex")

# mutant log2FoldChange bar plot with log2FC cutoff of 2, -2 & padjust < 0.05
sig_res <- sig_res[order(sig_res$log2FoldChange),]
ggplot(data = sig_res, aes(y = reorder(feature, sort(as.numeric(log2FoldChange))), x = log2FoldChange, fill = padj)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  labs(x = "Log2FoldChange (male/female)", y = "Pathways") +
  theme(
    text = element_text(size = 12), # Adjusts overall text size
    axis.title = element_text(size = 17), # Adjusts axis titles size
    axis.text = element_text(size = 16), # Adjusts axis text (tick labels) size
    legend.title = element_text(size = 14), # Adjusts legend title size
    legend.text = element_text(size = 14) # Adjusts legend text (labels) size
  )
# Mutant volcano plot (log2foldchange vs -log10 p-value) with log cutoff and p < 0.05

# Assuming sig_res is your dataframe and it contains log2FoldChange, pvalue, and description columns
# First, sort sig_res by pvalue to identify the most significant pathways
sig_res <- sig_res[order(sig_res$padj), ]

# You might want to label only a certain number of the most significant pathways, e.g., top 5
top_n <- 5
sig_res$label <- ifelse(seq_len(nrow(sig_res)) <= top_n, as.character(sig_res$feature), "")

ggplot(data = sig_res, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point() + 
  geom_label_repel(aes(label = label), box.padding = 0.35, point.padding = 0.2, 
                   size = 4, max.overlaps = Inf) +
  theme_minimal() +
  theme(
    text = element_text(size = 12), # Adjusts overall text size
    axis.title = element_text(size = 17), # Adjusts axis titles size
    axis.text = element_text(size = 16), # Adjusts axis text (tick labels) size
    legend.title = element_text(size = 14), # Adjusts legend title size
    legend.text = element_text(size = 14) # Adjusts legend text (labels) size
  )


### Heatmap with log2FC cutoff of 1.5, -1.5 for both genotypes (control vs. mutant)
source("DESeq2_function.R")

# Run the function on your own data
heat_res =  DEseq2_function(abundance_data_filtered, metadata, "Sex")
heat_res$feature =rownames(heat_res)
heat_res_desc = inner_join(heat_res,metacyc_daa_annotated_results_df, by = "feature")
heat_res_desc = heat_res_desc[, -c(8:13)]
View(heat_res_desc)

# Filter to only include significant pathways
heat_res_desc$padj <- as.numeric(heat_res_desc$padj)
heat_res_desc$log2FoldChange <- as.numeric(heat_res_desc$log2FoldChange)

heat_sig_res = heat_res_desc %>%
  filter(padj < 0.05 & (log2FoldChange > 1.5 | log2FoldChange < -1.5))

colnames(heat_sig_res)[7] <- "old_feature"
colnames(heat_sig_res)[8] <- "feature"

heat_abundance_desc <- abundance_desc[abundance_desc$feature %in% heat_sig_res$feature, ]

rownames(heat_abundance_desc) = NULL

# Generate a heatmap with log2FC cutoff of 1.5, -1.5 for BOTH genotypes
"feature" %in% colnames(heat_abundance_desc) #check
nrow(heat_abundance_desc) > 0 #check

# Combine 'Sex' and 'Genotype' into a new column 'Group'
sex_gen_metadata <- metadata %>% 
  mutate(Group = paste(Sex, Genotype, sep = "_"))

pathway_heatmap(abundance = heat_abundance_desc %>% column_to_rownames("feature"), 
                metadata = sex_gen_metadata, 
                group = "Group")

# ### pathway errorbar
# selected_pathways <- metacyc_daa_annotated_results_df %>% filter(p_adjust < 0.05) %>% select(c("feature","p_adjust"))
# reduced_metacyc_daa_annotated_results_df <- metacyc_daa_annotated_results_df[metacyc_daa_annotated_results_df$feature %in% selected_pathways$feature, ]
# 
# 
# p <- pathway_errorbar(abundance = abundance_data_filtered %>% column_to_rownames("pathway"),
#                       daa_results_df = reduced_metacyc_daa_annotated_results_df,
#                       Group = metadata$Sex,
#                       ko_to_kegg = FALSE,
#                       p_values_threshold = 0.05,
#                       order = "group",
#                       select = NULL,
#                       p_value_bar = TRUE,
#                       colors = NULL,
#                       x_lab = "feature")
# 
# # Reshape the data from wide to long format
# long_format <- pivot_longer(heat_abundance_desc, cols = -feature, names_to = "sample", values_to = "abundance")
# 
# # Convert abundance to percentage (if needed)
# # For example, if you want to convert raw abundance values to percentages of their row sums:
# long_format$abundance <- long_format$abundance / sum(long_format$abundance) * 100
# 
# # Create the boxplot
# ggplot(long_format, aes(x = feature, y = abundance)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x labels if needed
#   labs(x = "Pathway", y = "Abundance (%)", title = "Pathway Abundance Boxplots") +
#   scale_y_continuous(labels = percent_format()) # Format y-axis labels as percentages
# 
# boxplot(heat_abundance_desc$feature)
