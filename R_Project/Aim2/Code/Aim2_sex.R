# load libraries
library(tidyverse)
library(phyloseq)
library(ggsignif)
library(picante)
library(ggplot2)
library(RColorBrewer)


#load data and subset

load("FD_rare_only_separate.Rdata")

#duplicate phylosesq object
phyloseq_FD<- FD_rare_only_separate 

#extract the sample data
sample_phyloseq_FD <- sample_data(phyloseq_FD)

# Retain specific column
phyloseq_FD_subset <- sample_phyloseq_FD [, c("Sample.Name", "Age.New.Bin", "Cage.ID", "Experiment.Group", 
                                              "Genotype", "Mouse.ID", "Phenotype.score", 
                                              "FD.severity", "Sex", "Weight.grams")]

# Replace the sample data in your phyloseq object with the updated one
sample_data(FD_rare_only_separate) <- phyloseq_FD_subset

FD_phyloseq_subset <- FD_rare_only_separate

# extract alpha diversity information
alphadiv <- estimate_richness(FD_phyloseq_subset)
samp_dat <- sample_data(FD_phyloseq_subset)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)


#### RELATIVE ABUNDANCE ####
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
                                      measure.vars = 2:34, # melt using all sample columns
                                      variable.name = "Sample", # this will be the default x axis name
                                      value.name = "Relative_Abundance", # this will be the default y axis name
                                      as.is = TRUE # prevents type conversion of strings to factors
) #16929     2

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
  facet_nested(~ factor(Sex, , levels=c("Male", "Female")) + Genotype,scales = "free_x")+
  #facet_grid(~factor(Sex, , levels=c("Male", "Female")), scales = "free_x")+ #use this if you dont want to include genotype
  ylab('Relative abundance (%)')+
  xlab("Sample name")+
  facet_grid(~factor(Sex, , levels=c("Male", "Female")), scales = "free_x")+
  #  ggtitle("Relative abundance of phylum across mouse sexes")+
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(angle=+90, vjust=0.5, hjust=0))

ra.p



#### ALPHA ####

##### SHANNON #####

# run Wilxocon test on shannon indexes
wilcox_sh <- wilcox.test(Shannon ~ Sex, data=samp_dat_wdiv, exact = FALSE)
wilcox_sh

# graph shannon indexes to sex levels and annotate significance
shannon_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=Shannon)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  ggtitle("") +
  ylab('Shannon index') +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(4.4),
              annotations = c("0.5231"), textsize = 4) +
  theme_classic(base_size = 18)+
  scale_color_brewer(palette = "Set2")
shannon_sex

##### CHAO1 #####

# run Wilxocon test on Chao1 indexes
wilcox_ch <- wilcox.test(Chao1 ~ Sex, data=samp_dat_wdiv, exact = FALSE)
wilcox_ch

# graph Chao1 indexes to sex levels and annotate significance
chao1_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=Chao1)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  #  ggtitle("Chao1 and Sex") +
  ylab('Chao1 index') +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(180),
              annotations = c("0.5375"), textsize = 5) +
  theme_classic(base_size = 18)+
  scale_color_brewer(palette = "Set2")
chao1_sex


##### FAITH PHYLOGENETIC DIVERSITY #####
# import tree
tree <- read_tree("FD_tree.nwk")

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(FD_rare_only_separate)), phy_tree(tree), include.root=T) 

# add the PD column to subseted phyloseq object
sample_data(FD_rare_only_separate)$PD <- phylo_dist$PD

# extract the PD column from the sample data
PD_column <- sample_data(FD_rare_only_separate)$PD

# Add the PD column to the working data frame
samp_dat_wdiv$PD <- PD_column

# run Wilxocon test on pd indexes
wilcox_pd <- wilcox.test(PD ~ Sex, data=samp_dat_wdiv, exact = FALSE)
wilcox_pd

# graph pd indexes to sex levels and annotate significance
pd_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=PD)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  #  ggtitle("Faith's Phylogenetic Diversity and Sex") +
  ylab("Faith's Phylogenetic Diversity index") +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(15.3),
              annotations = c("0.3007"), textsize = 5) +
  theme_classic(base_size = 18)+
  scale_color_brewer(palette = "Set2")
pd_sex


##### PIELOU'S EVENNESS #####

# look for the total number of species in the subsetted data
total_species <- as.numeric(ntaxa(FD_rare_only_separate))

# calculate pielou's evenness index
samp_dat_wdiv$pielou <- samp_dat_wdiv$Shannon/log(total_species)

# run Wilxocon test on pielou's evenness indexes
wilcox_pielou <- wilcox.test(pielou ~ Sex, data=samp_dat_wdiv, exact = FALSE)
wilcox_pielou

# graph pd indexes to sex levels and annotate significance
pielou_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=pielou)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  #  ggtitle("Pielou's Evenness Index and Sex") +
  ylab("Pielou's Evenness index") +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(0.71),
              annotations = c("0.5231"), textsize = 5) +
  theme_classic(base_size = 18)+
  scale_color_brewer(palette = "Set2")
pielou_sex


#### BETA DIVERSITY ####

##### BRAY CURTIS DISTANCE #####
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

# perform permanova 
dist.bray <- vegdist(t(otu_table(FD_phyloseq_subset)), method="bray") # Bray-curtis

# view permanova results
sex_perm_bray <- adonis2(dist.bray ~ Sex,data = FD_phyloseq_sam_df)
sex_perm_bray



##### JACCARD DISTANCE #####

# create ordination for bray curtis
sex.jaccard.ord = ordinate(FD_phyloseq_subset, "PCoA", "jaccard", na.rm = TRUE) 

# plot ordination
sex.jaccard.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                    sex.jaccard.ord, #change this
                                    color = "Sex")
#                                    title = "Jaccard index across mouse sexes") #change this 

# format your plot
sex.jaccard.ord.pp = sex.jaccard.ord.p + #change this
  geom_point(size = 4)+ 
  scale_color_discrete(name = "Mouse sex")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
sex.jaccard.ord.pp

# perform permanova
dist.jaccard <- vegdist(t(otu_table(FD_phyloseq_subset)), method="jaccard")

# view permanova results
sex_perm_jaccard <- adonis2(dist.jaccard ~ Sex, data=FD_phyloseq_sam_df)




##### WEIGHTED UNIFRAC DISTANCE #####

#create ordination for bray curtis
sex.wunifrac.ord = ordinate(FD_phyloseq_subset, method = "PCoA", distance = "wunifrac") 

#plot ordination
sex.wunifrac.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                     sex.wunifrac.ord, #change this
                                     color = "Sex") #change this to the age category
#                                     title = "Weighted UniFrac distance across mouse sexes") #change this

#format your plot
sex.wunifrac.ord.pp= sex.wunifrac.ord.p + #change this
  geom_point(size = 4) +
  scale_color_discrete(name = "Mouse sex")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
sex.wunifrac.ord.pp

#perform permanova
dist.wunifrac <- UniFrac(FD_phyloseq_subset, weighted=TRUE) # Weighted UniFrac

# view permanova results
sex_perm_wunifrac <- adonis2(dist.wunifrac ~ Sex, data=FD_phyloseq_sam_df)




##### UNWEIGHTED UNIFRAC DISTANCE #####

# create ordination for unifrac
sex.unifrac.ord = ordinate(FD_phyloseq_subset, method = "PCoA", distance = "unifrac", weighted = FALSE)

# plot ordination
sex.unifrac.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                    sex.unifrac.ord, #change this
                                    color = "Sex") #change this to the age category
#                                    title = "Unweighted UniFrac distance across mouse sexes") #change this

# format your plot
sex.unifrac.ord.pp= sex.unifrac.ord.p + #change this
  geom_point(size = 4)+ 
  scale_color_discrete(name = "Mouse sex")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
sex.unifrac.ord.pp

# perform PERMANOVA
dist.unifrac <- UniFrac(FD_phyloseq_subset, weighted=FALSE) # Unweighted UniFrac

# view permanova results
sex_perm_unifrac <- adonis2(dist.unifrac ~ Sex, data=FD_phyloseq_sam_df)


#### EXPORT PLOTS AS PNGS ####

# Relative abundance one is at the end of part "RELATIVE ABUNDANCE"

png(("aim2_sex_relative_abundance.png"), width = 12, height = 5, units = "in", res = 300)
ra.p
dev.off()

png(("aim2_sex_shannon.png"), width = 5, height = 5, units = "in", res = 300)
shannon_sex
dev.off()

png(("aim2_sex_chao1.png"), width = 5, height = 5, units = "in", res = 300)
chao1_sex
dev.off()

png(("aim2_sex_pd.png"), width = 5, height = 5, units = "in", res = 300)
pd_sex
dev.off()

png(("aim2_sex_pielou.png"), width = 5, height = 5, units = "in", res = 300)
pielou_sex
dev.off()

png(("aim2_sex_bray.png"), width = 7.5, height = 5, units = "in", res = 300)
sex.bray.ord.pp
dev.off()

png(("aim2_sex_jaccard.png"), width = 7.5, height = 5, units = "in", res = 300)
sex.jaccard.ord.pp
dev.off()

png(("aim2_sex_wunifrac.png"), width = 7.5, height = 5, units = "in", res = 300)
sex.wunifrac.ord.pp
dev.off()

png(("aim2_sex_unifrac.png"), width = 7.5, height = 5, units = "in", res = 300)
sex.unifrac.ord.pp
dev.off()


#### SUMMARY: Wilcoxon results ####
wilcox_sh

wilcox_ch

wilcox_pd

wilcox_pielou

sex_perm_bray

sex_perm_jaccard

sex_perm_wunifrac

sex_perm_unifrac
