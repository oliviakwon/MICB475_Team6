##### Install packages #####
# Start by installing all necessary packages when asked if you want to install
# from source, please just type Yes in the terminal below

# If you don't have BiocManager, here is the code to install it
# A lot of you probably already have this so you can skip
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# when asked if you want to update all, some or none please type "n" for none

# After installing all of its above dependencies, install ggpicrust2
install.packages("ggpicrust2")
install.packages("readxl")

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

#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "picrust2_out_pipeline/pathabun_exported/pathway_abundance.tsv"
#Skipped first row (passage: exported from biome file)
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
abundance_data = as.data.frame(abundance_data)

#Import your metadata file, no need to filter yet
metadata <- read_excel("dysautonomia_metadata.xlsx")

#Looking at Age.months
#Since Age.months category has 3 variants
# filter the metadata to include only 2 at a time 
# Cross comparison: (3 months vs. 6 months), (3 months vs. 9 months), (6 months vs. 9 months)
# three_six_months_metadata <- metadata[metadata$Age.Months == c('3 months', '6 months'),]
# three_nine_months_metadata <- metadata[metadata$Age.Months == c('3 months', '9 months'),]
# six_nine_months_metadata <- metadata[metadata$Age.Months == c('6 months', '9 months'),]

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$Sex),]

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = metadata$'sample-id'
sample_names = append(sample_names, "#OTU ID")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering
# abundance_data_filtered[ , 'pathway'] = NA


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

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

colnames(abundance_data_filtered)[colnames(abundance_data_filtered) == "#OTU ID"] <- "pathway"

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
num_columns <- ncol(abundance_desc)
print(num_columns) #221 so 221 - 6("feature, "method", "group1", "group2", "p_values", "adj_method", "p_adjust", "description") = 213
abundance_desc = abundance_desc[,-c(215:ncol(abundance_desc))] 

# Generate a heatmap
"feature" %in% colnames(abundance_desc) #check
nrow(abundance_desc) > 0 #check

pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata, group = "Sex")

# Generate pathway PCA plot
pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadata, group = "Sex")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")

# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata, "Sex")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)
# filter by Log2fold change since we have a lot of pathways
# keep only greater than 2 fold change or less than -2
# sig_res = res_desc %>%
#   filter(log2FoldChange > 2 | log2FoldChange < -2)

sig_res <- sig_res[order(sig_res$log2FoldChange),]
ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")


# Volcano plot (log2foldchange vs -log10 p-value)

# Assuming sig_res is your dataframe and it contains log2FoldChange, pvalue, and description columns
# First, sort sig_res by pvalue to identify the most significant pathways
sig_res <- sig_res[order(sig_res$pvalue), ]

# You might want to label only a certain number of the most significant pathways, e.g., top 5
top_n <- 5
sig_res$label <- ifelse(seq_len(nrow(sig_res)) <= top_n, as.character(sig_res$description), "")

ggplot(data = sig_res, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point() + 
  geom_label_repel(aes(label = label), box.padding = 0.35, point.padding = 0.2, 
                   size = 3, max.overlaps = Inf) +
  theme_minimal()

