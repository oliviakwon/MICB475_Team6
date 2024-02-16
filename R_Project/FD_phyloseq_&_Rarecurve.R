#!/usr/bin/env Rscript
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

#### Load data ####
# Change file paths as necessary
metafp <- "QIIME2/export/FD_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "QIIME2/export/FD_feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "QIIME2/export/FD_taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "QIIME2/export/FD_tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df_meta <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df_meta)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP_META <- sample_data(samp_df_meta)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>% # remove confidence column
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
FD <- phyloseq(OTU, SAMP_META, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(FD)
sample_data(FD)
tax_table(FD)
phy_tree(FD)

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
FD_filt <- subset_taxa(FD,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
FD_filt_nolow <- filter_taxa(FD_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
FD_filt_nolow_samps <- prune_samples(sample_sums(FD_filt_nolow)>100, FD_filt_nolow)

FD_final <- (FD_filt_nolow_samps) 

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(FD_final))), cex=0.1)
FD_rare <- rarefy_even_depth(FD_final, rngseed = 1, sample.size = 9069)# sample size look at table.qzv (sampling depth)


##### Saving #####
save(FD_final, file="FD_final.RData")
save(FD_rare, file="FD_rare.RData")
