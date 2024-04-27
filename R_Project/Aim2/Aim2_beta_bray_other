# load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tibble)
library(vegan)
library(PMCMR)

# load data and subset

load("FD_rare_only_separate.Rdata")

#duplicate phylosesq object
phyloseq_FD<- FD_rare_only_separate 

#extract the sample data
sample_phyloseq_FD <- sample_data(phyloseq_FD)

# Retain specific column
phyloseq_FD_subset <- sample_phyloseq_FD [, c("Sample.Name", "Age.New.Bin", "Cage.ID", "Experiment.Group", 
                                              "Genotype", "Mouse.ID", "Phenotype.score", 
                                              "FD.severity", "Sex", "Weight.grams")]

#make sample data into a dataframe
phyloseq_FD_sam = sample_data(phyloseq_FD_subset)
phyloseq_FD_sam_df <- data.frame(phyloseq_FD_sam)

#set up quantile for mouse weight
phyloseq_FD_sam_df1 <- phyloseq_FD_sam_df %>%
  mutate(weight.quartile = ntile(Weight.grams, 3)) 

#define quantile for mouse weight
# Add a new column based on weight categories
phyloseq_FD_sam_df2 <- phyloseq_FD_sam_df1 %>%
  mutate(Weight_Category = ifelse(weight.quartile == 3, "High",
                                  ifelse(weight.quartile == 2, "Medium", "Low")))

as.factor(phyloseq_FD_sam_df2$Weight_Category)

# Replace the sample data in your phyloseq object with the updated one
sample_data(FD_rare_only_separate) <- phyloseq_FD_sam_df2

FD_phyloseq_subset <- FD_rare_only_separate

##### BETA DIVERSITY: BRAY CURTIS #####

###### Sex ######
#plot the order in Male, Female
sample_data(FD_phyloseq_subset)$Sex <- factor((sample_data(FD_phyloseq_subset)$Sex),
                                              levels = c("Male", "Female"))

#create ordination for bray curtis
sex.bray.ord = ordinate(FD_phyloseq_subset, "PCoA", "bray", na.rm = TRUE) 

#plot ordination
sex.bray.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                 sex.bray.ord, #change this
                                 color = "Sex") #change this to the age category
#                                 title = "Bray-Curtis dissimilarity across mouse sexes") #change this

#format your plot
sex.bray.ord.pp = sex.bray.ord.p + #change this
  geom_point(size = 4)+ 
  scale_color_discrete(name = "Mouse sex")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
sex.bray.ord.pp

# view permanova results
sex_perm_bray <- adonis2(dist.bray ~ Sex,data = FD_phyloseq_sam_df)
sex_perm_bray


###### Genotype ######

#plot the order in Male, Female
sample_data(FD_phyloseq_subset)$Genotype <- factor((sample_data(FD_phyloseq_subset)$Genotype),
                                                   levels = c("Control", "Mutant"))

#create ordination for bray curtis
genotype.bray.ord = ordinate(FD_phyloseq_subset, "PCoA", "bray", na.rm = TRUE) 

#plot ordination
genotype.bray.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                      genotype.bray.ord, #change this
                                      color = "Genotype") #change this to the age category
#                                      title = "Bray-Curtis dissimilarity across mouse genotypes") #change this

#format your plot
genotype.bray.ord.pp = genotype.bray.ord.p + #change this
  geom_point(size = 4)+ 
  scale_color_discrete(name = "Mouse genotype")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
genotype.bray.ord.pp

# perform permanova 
#dist.bray <- vegdist(t(otu_table(FD_phyloseq_subset)), method="bray") # Bray-curtis

#FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
#FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results
genotype_perm_bray <- adonis2(dist.bray ~ Genotype,data=FD_phyloseq_sam_df)
genotype_perm_bray


###### Weight ######

#plot the order in low, middle, high
sample_data(FD_phyloseq_subset)$Weight_Category <- factor((sample_data(FD_phyloseq_subset)$Weight_Category),
                                                          levels = c("High", "Medium", "Low"))

#create ordination for bray curtis
weight.bray.ord = ordinate(FD_phyloseq_subset, "PCoA", "bray", na.rm = TRUE) 

#plot ordination
weight.bray.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                    weight.bray.ord, #change this
                                    color = "Weight_Category") #change this to the weight category
#                                    title = "Bray-Curtis dissimilarity across mouse weight")

weight.cat.bray.ord.pp= weight.bray.ord.p + #change this
  geom_point(size = 4)+ 
  scale_color_discrete(name = "Mouse weight") +
  labs(color = "Weight category") +
  stat_ellipse(type = "norm") +
  theme_classic(base_size = 18) +
  scale_color_brewer(palette = "Set2")
weight.cat.bray.ord.pp

# perform permanova 
dist.bray <- vegdist(t(otu_table(FD_phyloseq_subset)), method="bray") # Bray-curtis

FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results
weight_cat_perm_bray <- adonis2(dist.bray ~ Weight_Category,data=FD_phyloseq_sam_df)
weight_cat_perm_bray



#### EXPORT PLOTS AS PNGS ####

png(("aim2_sex_bray.png"), width = 7.5, height = 5, units = "in", res = 300)
sex.bray.ord.pp
dev.off()

png(("aim2_genotype_bray.png"), width = 7.5, height = 5, units = "in", res = 300)
genotype.bray.ord.pp
dev.off()

png(("aim2_weight_categorical.bray.png"), width = 7.5, height = 5, units = "in", res = 300)
weight.cat.bray.ord.pp
dev.off()


#### SUMMARY: permanova and post-hoc results ####
sex_perm_bray

genotype_perm_bray

weight_perm_bray
