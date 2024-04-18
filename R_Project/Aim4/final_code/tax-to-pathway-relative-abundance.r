
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

# Step 1: Use tax_glom to calculate counts belonging to the same genus
sample_data(FD_final_only_separate) <- phyloseq_FD_subset
FD_phyloseq_subset <- FD_final_only_separate %>% tax_glom('Genus')

# Step 2: Merge OTU table and Taxonomy table from phyloseq_FD_subset 
taxa_abundances <- psmelt(FD_phyloseq_subset)
any(is.na(taxa_abundances$Genus)) # FALSE
# Replace any 'NA' in Genus with the Family value
# taxa_abundances$Genus <- ifelse(is.na(taxa_abundances$Genus), as.character(taxa_abundances$Family), taxa_abundances$Genus)

# Keep only mutant genotypes
taxa_abundances <- taxa_abundances[taxa_abundances$Genotype == 'Mutant', ]

#Importing the pathway abundance file from PICrsut2
abundance_file <- "pathway_abundance.tsv"
#Skipped first row (passage: exported from biome file)
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
abundance_data = as.data.frame(abundance_data)

sample_names = taxa_abundances$'Sample'
sample_names = append(sample_names, "#OTU ID")
abundance_data = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering


## relative abundance from abundance_data
# Convert to long format 
pathway_rel_abundance <- abundance_data %>%
  pivot_longer(
    cols = -"#OTU ID", 
    names_to = "Sample", 
    values_to = "Count" 
  )

# Calculate the total counts for each sample
totals <- pathway_rel_abundance %>%
  group_by(Sample) %>%
  summarise(TotalCount = sum(Count))

# Join the total count back to the original data frame
pathway_rel_abundance <- pathway_rel_abundance %>%
  left_join(totals, by = "Sample")

# Calculate relative abundance
pathway_rel_abundance <- pathway_rel_abundance %>%
  mutate(RelativeAbundance = Count / TotalCount)

# Remove the TotalCount column, pathway_rel_abundance contains 
# pathway relative abundance information for each sample in each pathway
pathway_rel_abundance <- pathway_rel_abundance %>%
  select(-TotalCount)

## Step 3: Sum up the relative abundances for each Genus
# Calculate the total counts for each sample to use for relative abundance
total_counts <- taxa_abundances %>% 
  group_by(Sample) %>% 
  summarise(Total = sum(Abundance))

# Join total counts back to the taxa abundance data
taxa_abundances <- taxa_abundances %>%
  left_join(total_counts, by = "Sample")

# Calculate relative abundance
taxa_abundances <- taxa_abundances %>%
  mutate(RelativeAbundance = Abundance / Total)

# Relative abundance for each genus in each sample
relative_abundance_df <- taxa_abundances %>%
  group_by(Sample, Genus) %>%
  summarise(RelativeAbundance = sum(RelativeAbundance), .groups = 'drop') %>%
  spread(key = Genus, value = RelativeAbundance)

# Relative abundance for each genus in each sample in long format
relative_abundance_long_df <- relative_abundance_df %>%
  gather(key = "Genus", value = "RelativeAbundance", -Sample)

#Rename columns to avoid confusion between relative abundances
colnames(pathway_rel_abundance)[colnames(pathway_rel_abundance) == "RelativeAbundance"] <- "PathwayRelativeAbundance"
colnames(pathway_rel_abundance)[colnames(pathway_rel_abundance) == "#OTU ID"] <- "Pathways"
colnames(relative_abundance_long_df)[colnames(relative_abundance_long_df) == "RelativeAbundance"] <- "GenusRelativeAbundance"

# Pathway abundance into categorical values for the plot
# Using K-means Clustering, determine breaks based on quantiles
quantile_breaks <- quantile(pathway_rel_abundance$PathwayRelativeAbundance, probs = seq(0, 1, by = 0.25), na.rm = TRUE)

category_labels <- c("1st Quartile", "2nd Quartile", "3rd Quartile", "4th Quartile")

# Use the labels in the cut function to make sure that the resulting variable is a factor
pathway_rel_abundance$AbundanceCategory <- cut(
  pathway_rel_abundance$PathwayRelativeAbundance,
  breaks = quantile_breaks,
  include.lowest = TRUE,
  labels = category_labels
)

# The 4 pathways with -2 < log2FC > 2, padjust < 0.05
pathway_rel_abundance <- pathway_rel_abundance %>%
  filter(Pathways %in% c("PWY-4722", "PWY-5419", "PWY-5420", "PWY-7007"))

library(dplyr)

pathway_rel_abundance <- pathway_rel_abundance %>%
  mutate(feature = case_when(
    Pathways == "PWY-4722" ~ "creatinine degradation II",
    Pathways == "PWY-5419" ~ "catechol degradation to 2-oxopent-4-enoate II",
    Pathways == "PWY-5420" ~ "catechol degradation II (meta-cleavage pathway)",
    Pathways == "PWY-7007" ~ "methyl ketone biosynthesis",
    TRUE ~ NA_character_  # if no match
  ))


# Merge the dataframes by 'Sample'
merged_data <- merge(relative_abundance_long_df, pathway_rel_abundance, by = "Sample")

# Pathway = x axis, Genus is = y axis
ggplot(merged_data, aes(x = feature, y = Genus, size = AbundanceCategory, color = GenusRelativeAbundance)) +
  geom_point(alpha = 0.7) +
  scale_size_manual(values = c(3, 6, 9, 12)) +  # sizes for categories
  scale_color_gradient(low = "purple", high = "yellow") +
  theme_minimal() +
  labs(x = "Pathways", y = "Genus", color = "Genus Relative Abundance", size = "Pathway Abundance") +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),    # Legend titles
        legend.text = element_text(size = 20),     # Legend text
        axis.title = element_text(size = 20),      # Axis titles
        axis.text.y = element_text(size = 14),     # y-axis labels
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1))  # x-axis labels

