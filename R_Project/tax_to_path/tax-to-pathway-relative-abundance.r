
#---------- load libraries----------
library(tidyverse)
library(phyloseq)
library(ggsignif)
library(picante)
library(ggh4x)


#----------load data and subset----------

load("FD_final_only_separate.Rdata")

#duplicate phylosesq object
phyloseq_FD<- FD_final_only_separate 

#extract the sample data
sample_phyloseq_FD <- sample_data(phyloseq_FD)

# Retain specific column
phyloseq_FD_subset <- sample_phyloseq_FD [, c("Sample.Name", "Age.New.Bin", "Cage.ID", "Experiment.Group", 
                                              "Genotype", "Mouse.ID", "Phenotype.score", 
                                              "FD.severity", "Sex", "Weight.grams")]

# Count the number of rows where the value in the "sex" column is "female"
num_females <- sum(phyloseq_FD_subset$Sex == 'Female') #43 female

num_females_control <- sum(phyloseq_FD_subset$Sex == 'Female' &  #22
                             phyloseq_FD_subset$Genotype == 'Control')

num_females_mutant <- sum(phyloseq_FD_subset$Sex == 'Female' &  #21
                            phyloseq_FD_subset$Genotype == 'Mutant')

num_males <- sum(phyloseq_FD_subset$Sex == 'Male') #16 female

num_Male_control <- sum(phyloseq_FD_subset$Sex == 'Male' &  #3
                          phyloseq_FD_subset$Genotype == 'Control')

num_male_mutant <- sum(phyloseq_FD_subset$Sex == 'Male' &  #13
                         phyloseq_FD_subset$Genotype == 'Mutant')



# Replace the sample data in your phyloseq object with the updated one
sample_data(FD_final_only_separate) <- phyloseq_FD_subset

FD_phyloseq_subset <- FD_final_only_separate

#---------- RELATIVE ABUNDANCE ----------
#only need otu table as dataframe
FD_phyloseq_otu<- otu_table(FD_phyloseq_subset)
FD_otu<-as.data.frame(FD_phyloseq_otu)

#convert abs read to relative abundance 
FD_otu_relative <- FD_otu %>%
  dplyr::mutate_at(c(1:ncol(FD_otu)), funs(./sum(.)*100))

#test whether each column (aka sample) is 100%
sum(FD_otu_relative$SRR17661697.fastq.gz)

#change the row to column
FD_otu_relative1<- rownames_to_column(FD_otu_relative, var ="ASV")

#melt the dataframe
FD_otu_relative_melt<- reshape2::melt(data = FD_otu_relative1, 
                                      measure.vars = 2:60, # melt using all sample columns
                                      variable.name = "Sample", # this will be the default x axis name
                                      value.name = "Relative_Abundance", # this will be the default y axis name
                                      as.is = TRUE # prevents type conversion of strings to factors
) #65195     3

#combine taxa info to the FD_otu_relative_melt
FD_phyloseq_tax<- tax_table(FD_phyloseq_subset)
FD_tax<-as.data.frame(FD_phyloseq_tax)
FD_tax1<- rownames_to_column(FD_tax, var ="ASV")

test<- dplyr::left_join(FD_otu_relative_melt,FD_tax1, by = "ASV")

#subset sample data with sample name and age
FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

FD_phyloseq_sam_sub<- subset(FD_phyloseq_sam_df, select = c("Sample.Name", "Sex", "Genotype"))

FD_phyloseq_sam_sub1 <- rownames_to_column(FD_phyloseq_sam_sub, var = "Sample")

#combine test and FD_phyloseq_sam_sub1
taxa_ra<- dplyr::left_join(test,FD_phyloseq_sam_sub1, by = "Sample")

#clean up the phylum name
taxa_ra$Phylum<- gsub("^[a-z]__", "",taxa_ra$Phylum )

ra.p<- ggplot(data =taxa_ra, aes(x= Sample.Name, y = Relative_Abundance))+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  ylab('Relative abundance (%)')+
  xlab("Sample name")+
  facet_nested(~ factor(Sex, , levels=c("Male", "Female")) + Genotype,scales = "free_x")+
  #facet_grid(~factor(Sex, , levels=c("Male", "Female")), scales = "free_x")+
  ggtitle("Relative abundance of phylum across mouse sexes")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(angle=+90, vjust=0.5, hjust=0))

ra.p

#### Import files and preparing tables ####
#Importing the pathway abundance file from PICrsut2
abundance_file <- "pathway_abundance.tsv"
#Skipped first row (passage: exported from biome file)
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
abundance_data = as.data.frame(abundance_data)

## relative abundance from abundance_data
# Convert to long format 
pathway_rel_abundance <- abundance_data %>%
  pivot_longer(
    cols = -"#OTU ID", 
    names_to = "Sample", 
    values_to = "Count" 
  )

# Step 1: Use pivot_wider to transform the data frame such that 'Sample' becomes columns
taxa_ra_wide <- taxa_ra %>%
  pivot_wider(names_from = Sample, values_from = Relative_Abundance, 
              values_fill = list(Relative_Abundance = 0))

# Step 2: Replace 'NA' in Genus with the Family value
taxa_ra_wide$Genus <- ifelse(is.na(taxa_ra_wide$Genus), as.character(taxa_ra_wide$Family), taxa_ra_wide$Genus)

# # Step 3: Sum up the relative abundances for each Genus
# numeric_columns <- names(taxa_ra_wide)[12:ncol(taxa_ra_wide)]
# 
# genus_sums <- taxa_ra_wide %>%
#   group_by(Genus) %>%
#   summarise(across(all_of(numeric_columns), sum, na.rm = TRUE), .groups = "drop")
# 
# # Each row with the same genus has the total summed abundance for that genus
# taxa_ra_wide <- taxa_ra_wide %>%
#   left_join(genus_sums, by = "Genus")

# Step 4: Make sample names into a column called 'Sample'
taxa_ra_long <- taxa_ra_wide %>%
  pivot_longer(cols = names(taxa_ra_wide)[12:ncol(taxa_ra_wide)], names_to = "Sample", values_to = "Relative_Abundance")

# Step 5: Join with pathway relative abundance dataframe 
joined_df <- left_join(pathway_rel_abundance, taxa_ra_long, by = "Sample")

