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
experiment_data <- merged_data %>% group_by(Genus, Pathways) %>%
  group_modify(~broom::tidy(cor.test(.$GenusRelativeAbundance, .$PathwayRelativeAbundance, data = ., method = "spearman"))) %>%
  ungroup() %>%
  mutate(padjust = p.adjust(p.value))

padjust_0.05_experiment_data <- experiment_data %>% filter(padjust < 0.05)

# Scatter plot: filter merge data to filter only for species in padjust_0.05_experiment_data
# show padjust and estimate 

filtered_merged_data <- merged_data[merged_data$Genus == "g__Romboutsia", ]
colnames(filtered_merged_data)[colnames(filtered_merged_data) == "Pathways"] <- "old_pathways"
colnames(filtered_merged_data)[colnames(filtered_merged_data) == "feature"] <- "Pathways"

filtered_merged_data <- filtered_merged_data %>%
  mutate(feature = case_when(
    Pathways == "PWY-4722" ~ "creatinine degradation II",
    Pathways == "PWY-5419" ~ "catechol degradation to 2-oxopent-4-enoate II",
    Pathways == "PWY-5420" ~ "catechol degradation II (meta-cleavage pathway)",
    Pathways == "PWY-7007" ~ "methyl ketone biosynthesis",
    TRUE ~ NA_character_  # if no match
  ))

ggplot(filtered_merged_data, aes(x = GenusRelativeAbundance, y = PathwayRelativeAbundance, color = Pathways)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm") +
  scale_color_brewer(palette = "Set2") +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", size = 6.5) +
  labs(x = "Genus Relative Abundance", y = "Pathway Relative Abundance") + 
  theme_classic() + 
  theme(legend.position = "right",
        legend.title = element_text(size = 18),    # Legend titles
        legend.text = element_text(size = 18),     # Legend text
        axis.title = element_text(size = 18),      # Axis titles
        axis.text.y = element_text(size = 18),     # y-axis labels
        axis.text.x = element_text(size = 18)) +    # x-axis labels 
  guides(color = guide_legend(nrow = 4))
