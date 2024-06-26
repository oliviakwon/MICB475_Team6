
#### Loading packages #### 
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(ape) 
library(vegan)

# finding indicator species for the significant variables 
# need to remove control samples? bc not FD affected 

#### Load physoleq object ####
load('FD_final_only_separate.RData') 

# rename the phyloseq object
meta_FD <- FD_final_only_separate
meta_FD_no_controls <- subset_samples(meta_FD, Genotype != "Control") 

#### Variable: Sex ####
## ISA Phylum ##
# Glomming the phyloseq object
meta_phylum <- tax_glom(meta_FD_no_controls, "Phylum", NArm = FALSE)
# convert counts into relative abundance
meta_phylum_RA <- transform_sample_counts(meta_phylum, fun=function(x) x/sum(x))
# calculate indicator values for all ASVs
isa_sex_phylum <- multipatt(t(otu_table(meta_phylum_RA)), cluster = sample_data(meta_phylum_RA)$'Sex')
summary(isa_sex_phylum)
# combine the taxanomy table with indicator species and indicator p values
taxtable <- tax_table(meta_FD_no_controls) %>% as.data.frame() %>% rownames_to_column(var="ASV")

sex_phylum <- isa_sex_phylum$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

## ISA Class ##
meta_class <- tax_glom(meta_FD_no_controls, "Class", NArm = FALSE)
meta_class_RA <- transform_sample_counts(meta_class, fun=function(x) x/sum(x))
isa_sex_class <- multipatt(t(otu_table(meta_class_RA)), cluster = sample_data(meta_class_RA)$'Sex')
summary(isa_sex_class)
taxtable <- tax_table(meta_FD_no_controls) %>% as.data.frame() %>% rownames_to_column(var="ASV")

sex_class <- isa_sex_class$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

## ISA Order ##
meta_order <- tax_glom(meta_FD_no_controls, "Order", NArm = FALSE)
meta_order_RA <- transform_sample_counts(meta_order, fun=function(x) x/sum(x))
isa_sex_order <- multipatt(t(otu_table(meta_order_RA)), cluster = sample_data(meta_order_RA)$'Sex')
summary(isa_sex_order)
taxtable <- tax_table(meta_FD_no_controls) %>% as.data.frame() %>% rownames_to_column(var="ASV")

sex_order <- isa_sex_order$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

## ISA Family ##
meta_family <- tax_glom(meta_FD_no_controls, "Family", NArm = FALSE)
meta_family_RA <- transform_sample_counts(meta_family, fun=function(x) x/sum(x))
isa_sex_family <- multipatt(t(otu_table(meta_family_RA)), cluster = sample_data(meta_family_RA)$'Sex')
summary(isa_sex_family)
taxtable <- tax_table(meta_FD_no_controls) %>% as.data.frame() %>% rownames_to_column(var="ASV")

sex_family <- isa_sex_family$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

## ISA Genus ##
meta_genus <- tax_glom(meta_FD_no_controls, "Genus", NArm = FALSE)
meta_genus_RA <- transform_sample_counts(meta_genus, fun=function(x) x/sum(x))
isa_sex_genus <- multipatt(t(otu_table(meta_genus_RA)), cluster = sample_data(meta_genus_RA)$'Sex')
summary(isa_sex_genus)
taxtable <- tax_table(meta_FD_no_controls) %>% as.data.frame() %>% rownames_to_column(var="ASV")

sex_genus <- isa_sex_genus$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

## ISA Species ##
meta_species <- tax_glom(meta_FD_no_controls, "Species", NArm = FALSE)
meta_species_RA <- transform_sample_counts(meta_species, fun=function(x) x/sum(x))
isa_sex_species <- multipatt(t(otu_table(meta_species_RA)), cluster = sample_data(meta_species_RA)$'Sex')
summary(isa_sex_species)
taxtable <- tax_table(meta_FD_no_controls) %>% as.data.frame() %>% rownames_to_column(var="ASV")

ISA_Sex_species <- isa_sex_species$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05, stat>0.80) %>% as.data.frame()
# no filter on stats = 14 indicator species
# stat> 0.80 = 8 indicator species
write.csv(ISA_Sex_species, "ISA_Sex_Species.csv", row.names = FALSE)

## ISA ASV ##
meta_ASV_RA <- transform_sample_counts(meta_FD_no_controls, fun=function(x) x/sum(x))
isa_sex_ASV <- multipatt(t(otu_table(meta_ASV_RA)), cluster = sample_data(meta_ASV_RA)$'Sex')
summary(isa_sex_ASV, indvalcomp=TRUE) 
taxtable <- tax_table(meta_FD_no_controls) %>% as.data.frame() %>% rownames_to_column(var="ASV")

ISA_Sex_ASV <- isa_sex_ASV$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05, stat>0.80) %>% as.data.frame()
# stat>0.85 yields 14 observed ASVs 
# stat>0.80 yields 21 observed ASVs

write.csv(ISA_Sex_ASV, "ISA_Sex_ASV.csv", row.names = FALSE)



