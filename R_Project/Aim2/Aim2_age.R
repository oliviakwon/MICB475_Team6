# load libraries
library(tidyverse)
library(phyloseq)
library(ggsignif)
library(picante)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tibble)
library(vegan)
library(PMCMR)

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


##### RELATIVE ABUNDANCE #####
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

FD_phyloseq_sam_sub<- FD_phyloseq_sam_df[, 1:2]

FD_phyloseq_sam_sub1 <- rownames_to_column(FD_phyloseq_sam_sub, var = "Sample")

#combine test and FD_phyloseq_sam_sub1
taxa_ra<- dplyr::left_join(test,FD_phyloseq_sam_sub1, by = "Sample")


ra.p<- ggplot(data =taxa_ra, aes(x= Sample.Name, y = Relative_Abundance))+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  ylab('Relative abundance (%)')+
  facet_grid(~factor(Age.New.Bin, levels=c("young", "middle", "old")), scales = "free_x")+
  ggtitle("Relative abundance of phylum across mouse age")+
  theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0))+
  theme_bw()
theme_update(strip.text.x = element_text(size = 28),
             strip.text.y = element_text(size = 28),
             axis.text.x = element_text(size = 24), # change to 18 for smaller x axis label
             axis.text.y = element_text(size = 24),
             axis.title.x = element_text(size = 28),
             axis.title.y = element_text(size = 28),
             legend.title = element_text(size = 28),
             legend.text = element_text(size = 28),
             plot.title = element_text(size = 28))


png(("Relative abundance of phylum across mouse age.png"), width = 1920, height = 1080)
ra.p
dev.off()


##### ALPHA DIVERSITY #####

###### SHANNON ######

# run Kruskal-wallis test on shannon indexes
kruskal_sh <- kruskal.test(Shannon ~ `Age.New.Bin`, data = samp_dat_wdiv)
kruskal_sh

# log the data and run ANOVA to see significance
lm_sh_vs_age_log <- lm(log(Shannon) ~ `Age.New.Bin`, data=samp_dat_wdiv)
anova_sh_vs_age_log <- aov(lm_sh_vs_age_log)
#summary(anova_sh_vs_age_log)
TukeyHSD(anova_sh_vs_age_log)

# graph shannon indexes to age levels and annotate significance
# based on the logged ANOVA above
shannon_a <- ggplot(samp_dat_wdiv, aes(x=`Age.New.Bin`, y=Shannon)) +
  scale_x_discrete(limits =c("young", "middle", "old")) +
  geom_boxplot() +
  geom_point(size = 5) +
  ggtitle("Shannon and Age") +
  ylab('Shannon index') +
  xlab('Mouse age') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("old","middle"), c("young", "middle"), c("young","old")),
              y_position = c(4.4, 4.4, 4.5),
              annotations = c("0.984","0.470","0.580"), textsize = 10) +
  theme_update(strip.text.x = element_text(size = 28),
               strip.text.y = element_text(size = 28),
               axis.text.x = element_text(size = 24), # change to 18 for smaller x axis label
               axis.text.y = element_text(size = 24),
               axis.title.x = element_text(size = 28),
               axis.title.y = element_text(size = 28),
               legend.title = element_text(size = 28),
               legend.text = element_text(size = 28),
               plot.title = element_text(size = 28))


###### CHAO1 ######

# run Kruskal-wallis test on Chao1 indexes
kruskal_ch <- kruskal.test(Chao1 ~ `Age.New.Bin`, data = samp_dat_wdiv)
kruskal_ch

# log the data and run ANOVA to see significance
lm_ch_vs_age_log <- lm(log(Chao1) ~ `Age.New.Bin`, data=samp_dat_wdiv)
anova_ch_vs_age_log <- aov(lm_ch_vs_age_log)
#summary(anova_ch_vs_age_log)
TukeyHSD(anova_ch_vs_age_log)

# graph Chao1 indexes to age levels and annotate significance
# based on the logged ANOVA above
chao1_a <- ggplot(samp_dat_wdiv, aes(x=`Age.New.Bin`, y=Chao1)) +
  scale_x_discrete(limits =c("young", "middle", "old")) +
  geom_boxplot() +
  geom_point(size = 5) +
  ggtitle("Chao1 and Age") +
  ylab('Chao1 index') +
  xlab('Mouse age') +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("old","middle"), c("young", "middle"), c("young","old")),
              y_position = c(175, 175, 185),
              annotations = c("0.883","0.377","0.901"), textsize = 10) +
  theme_update(strip.text.x = element_text(size = 28),
               strip.text.y = element_text(size = 28),
               axis.text.x = element_text(size = 24), # change to 18 for smaller x axis label
               axis.text.y = element_text(size = 24),
               axis.title.x = element_text(size = 28),
               axis.title.y = element_text(size = 28),
               legend.title = element_text(size = 28),
               legend.text = element_text(size = 28),
               plot.title = element_text(size = 28))


###### FAITH PHYLOGENETIC DIVERSITY ######
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

# run Kruskal-wallis test on Faith's PD
kruskal_pd <- kruskal.test(PD ~ `Age.New.Bin`, data = samp_dat_wdiv)
kruskal_pd

# log the data and run ANOVA to see significance
lm_pd_vs_age_log <- lm(log(PD) ~ `Age.New.Bin`, data=samp_dat_wdiv)
anova_pd_vs_age_log <- aov(lm_pd_vs_age_log)
#summary(anova_pd_vs_age_log)
TukeyHSD(anova_pd_vs_age_log)

# graph PD to age levels and annotate significance
# based on the logged ANOVA above
pd_a <- ggplot(samp_dat_wdiv, aes(x=`Age.New.Bin`, y=PD)) +
  scale_x_discrete(limits =c("young", "middle", "old")) +
  geom_boxplot() +
  geom_point(size = 5) +
  ggtitle("Faith's Phylogenetic Diversity and Age") +
  ylab('Phylogenetic diversity index') +
  xlab('Mouse age') +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("old","middle"), c("young", "middle"), c("young","old")),
              y_position = c(15.2, 15.2, 15.7),
              annotations = c("0.767","0.509","0.999"), textsize = 10) +
  theme_update(strip.text.x = element_text(size = 28),
               strip.text.y = element_text(size = 28),
               axis.text.x = element_text(size = 24), # change to 18 for smaller x axis label
               axis.text.y = element_text(size = 24),
               axis.title.x = element_text(size = 28),
               axis.title.y = element_text(size = 28),
               legend.title = element_text(size = 28),
               legend.text = element_text(size = 28),
               plot.title = element_text(size = 28))


###### PIELOU'S EVENNESS ######

# look for the total number of species in the subsetted data
total_species <- as.numeric(ntaxa(FD_rare_only_separate))

# calculate pielou's evenness index
samp_dat_wdiv$pielou <- samp_dat_wdiv$Shannon/log(total_species)

# run Kruskal-wallis test on Pielou's evenness indexes
kruskal_pielou <- kruskal.test(pielou ~ `Age.New.Bin`, data = samp_dat_wdiv)
kruskal_pielou

# log the data and run ANOVA to see significance
lm_pielou_vs_age_log <- lm(log(pielou) ~ `Age.New.Bin`, data=samp_dat_wdiv)
anova_pielou_vs_age_log <- aov(lm_pielou_vs_age_log)
#summary(anova_pielou_vs_age_log)
TukeyHSD(anova_pielou_vs_age_log)

# graph Pielou's evenness indexes to age levels and annotate significance
# based on the logged ANOVA above
pielou_a <- ggplot(samp_dat_wdiv, aes(x=`Age.New.Bin`, y=pielou)) +
  scale_x_discrete(limits =c("young", "middle", "old")) +
  geom_boxplot() +
  geom_point(size = 5) +
  ggtitle("Pielou's evenness index and age") +
  ylab("Pielou's evenness index") +
  xlab('Mouse age') +
  theme(axis.text.x=element_text(angle=0, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("old","middle"), c("young", "middle"), c("young","old")),
              y_position = c(0.71, 0.71, 0.73),
              annotations = c("0.984","0.470","0.580"), textsize = 10) +
  theme_update(strip.text.x = element_text(size = 28),
               strip.text.y = element_text(size = 28),
               axis.text.x = element_text(size = 24), # change to 18 for smaller x axis label
               axis.text.y = element_text(size = 24),
               axis.title.x = element_text(size = 28),
               axis.title.y = element_text(size = 28),
               legend.title = element_text(size = 28),
               legend.text = element_text(size = 28),
               plot.title = element_text(size = 28))



##### BETA DIVERSITY #####
# Capitalize levels in Age.New.Bin
sample_data(FD_phyloseq_subset)$Age.New.Bin <- str_to_title(sample_data(FD_phyloseq_subset)$Age.New.Bin)

# Print the updated column
print(sample_data(FD_phyloseq_subset)$Age.New.Bin)

###### BRAY CURTIS DISTANCE ######

#plot the order in young, middle, then old
sample_data(FD_phyloseq_subset)$Age.New.Bin <- factor((sample_data(FD_phyloseq_subset)$Age.New.Bin),
                                                      levels = c("Old", "Middle", "Young"))

#create ordination for bray curtis
age.bray.ord = ordinate(FD_phyloseq_subset, "PCoA", "bray", na.rm = TRUE) 

#plot ordination
age.bray.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                 age.bray.ord, #change this
                                 color = "Age.New.Bin") #change this to the age category
#                                 title = "Bray-Curtis dissimilarity across mouse ages" #change this

#format your plot
age.bray.ord.pp = age.bray.ord.p + #change this
  labs(color = "Age") +
  geom_point(size = 4) + 
  scale_color_discrete(name = "Mouse age")+
  theme_classic(base_size = 18) + 
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
age.bray.ord.pp

# perform permanova 
dist.bray <- vegdist(t(otu_table(FD_phyloseq_subset)), method="bray") # Bray-curtis

FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results
perm_bray <- adonis2(dist.bray ~ Age.New.Bin,data=FD_phyloseq_sam_df)

# post-hoc test for permanova
posthoc_bray <- RVAideMemoire::pairwise.perm.manova(resp = dist.bray, fact = sample_data(FD_phyloseq_subset)$Age.New.Bin, 
                                                    test = "Pillai", nperm = 999, progress = TRUE, p.method = "none")
df_bray <- reshape2::melt(posthoc_bray$p.value)
colnames(df_bray) <- c("Age1", "Age2", "pvalue")
df_bray <- df_bray[-which(is.na(df_bray$pvalue)), ]
df_bray$pvalue.adj <- p.adjust(p = df_bray$pvalue, method = "bonferroni", n = length(df_bray$pvalue))

#view post-hoc results
head(df_bray)



###### JACCARD DISTANCE ######

# create ordination for bray curtis
age.jaccard.ord = ordinate(FD_phyloseq_subset, "PCoA", 
                           "jaccard", na.rm = TRUE) 

# plot ordination
age.jaccard.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                    age.jaccard.ord, #change this
                                    color = "Age.New.Bin")
#                                    title = "Jaccard index across mouse ages" #change this

# format your plot
age.jaccard.ord.pp = age.jaccard.ord.p + #change this
  geom_point(size = 4) + 
  labs(color = "Age") +
  scale_color_discrete(name = "Mouse age")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm")  +
  scale_color_brewer(palette = "Set2")
age.jaccard.ord.pp

#perform permanova
dist.jaccard <- vegdist(t(otu_table(FD_phyloseq_subset)), method="jaccard")

FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results
perm_jaccard <- adonis2(dist.jaccard ~ Age.New.Bin, data=FD_phyloseq_sam_df)

# Post-hoc test for permanova
posthoc_jaccard <- RVAideMemoire::pairwise.perm.manova(resp = dist.jaccard, fact = sample_data(FD_phyloseq_subset)$Age.New.Bin, 
                                                       test = "Pillai", nperm = 999, progress = TRUE, p.method = "none")
df_jaccard <- reshape2::melt(posthoc_jaccard$p.value)
colnames(df_jaccard) <- c("Age1", "Age2", "pvalue")
df_jaccard <- df_jaccard[-which(is.na(df_jaccard$pvalue)), ]
df_jaccard$pvalue.adj <- p.adjust(p = df_jaccard$pvalue, method = "bonferroni", n = length(df_jaccard$pvalue))

# view post-hoc results
head(df_jaccard)


###### WEIGHTED UNIFRAC DISTANCE - AGE ######

#create ordination for bray curtis
age.wunifrac.ord = ordinate(FD_phyloseq_subset, method = "PCoA", distance = "wunifrac") 

#plot ordination
age.wunifrac.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                     age.wunifrac.ord, #change this
                                     color = "Age.New.Bin") #change this to the age category
#                                     title = "Weighted UniFrac distance across mouse ages") #change this

#format your plot
age.wunifrac.ord.pp= age.wunifrac.ord.p + #change this
  geom_point(size = 4) +
  labs(color = "Age") +
  scale_color_discrete(name = "Mouse age")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
age.wunifrac.ord.pp

#perform permanova
dist.wunifrac <- UniFrac(FD_phyloseq_subset, weighted=TRUE) # Weighted UniFrac

FD_phyloseq_sam <- sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

perm_wunifrac <- adonis2(dist.wunifrac ~ Age.New.Bin, data=FD_phyloseq_sam_df)

# Post-hoc test for permanova
posthoc_wunifrac <- RVAideMemoire::pairwise.perm.manova(resp = dist.wunifrac, fact = sample_data(FD_phyloseq_subset)$Age.New.Bin, 
                                                        test = "Pillai", nperm = 999, progress = TRUE, p.method = "none")
df_wunifrac <- reshape2::melt(posthoc_wunifrac$p.value)
colnames(df_wunifrac) <- c("Age1", "Age2", "pvalue")
df_wunifrac <- df_wunifrac[-which(is.na(df_wunifrac$pvalue)), ]
df_wunifrac$pvalue.adj <- p.adjust(p = df_wunifrac$pvalue, method = "bonferroni", n = length(df_wunifrac$pvalue))

# view post-hoc results
head(df_wunifrac)



###### UNWEIGHTED UNIFRAC DISTANCE - AGE ######

# create ordination for unifrac
age.unifrac.ord = ordinate(FD_phyloseq_subset, method = "PCoA", distance = "unifrac", weighted = FALSE)

# plot ordination
age.unifrac.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                    age.unifrac.ord, #change this
                                    color = "Age.New.Bin") #change this to the age category
#                                    title = "Unweighted UniFrac distance across mouse ages") #change this

# format your plot
age.unifrac.ord.pp= age.unifrac.ord.p + #change this
  geom_point(size = 4) +
  labs(color = "Age") +
  scale_color_discrete(name = "Mouse age")+
  theme_classic(base_size = 18) +
  stat_ellipse(type = "norm") +
  scale_color_brewer(palette = "Set2")
age.unifrac.ord.pp

# perform PERMANOVA
dist.unifrac <- UniFrac(FD_phyloseq_subset, weighted=FALSE) # Unweighted UniFrac

FD_phyloseq_sam <- sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results
perm_unifrac <- adonis2(dist.unifrac ~ Age.New.Bin, data=FD_phyloseq_sam_df)

# post-hoc test for permanova
posthoc_unifrac <- RVAideMemoire::pairwise.perm.manova(resp = dist.unifrac, fact = sample_data(FD_phyloseq_subset)$Age.New.Bin, 
                                                       test = "Pillai", nperm = 999, progress = TRUE, p.method = "none")
df_unifrac <- reshape2::melt(posthoc_unifrac$p.value)
colnames(df_unifrac) <- c("Age1", "Age2", "pvalue")
df_unifrac <- df_unifrac[-which(is.na(df_unifrac$pvalue)), ]
df_unifrac$pvalue.adj <- p.adjust(p = df_unifrac$pvalue, method = "bonferroni", n = length(df_unifrac$pvalue))

# view post-hoc results
head(df_unifrac)

##### EXPORT PLOTS AS PNGS #####

png(("aim2_age_relative_abundance.png"), width = 1000, height = 800)
ra.p
dev.off()

png(("aim2_age_shannon.png"), width = 1000, height = 800)
shannon_a
dev.off()

png(("aim2_age_chao1.png"), width = 1000, height = 800)
chao1_a
dev.off()

png(("aim2_age_pd.png"), width = 1000, height = 800)
pd_a
dev.off()

png(("aim2_age_pielou.png"), width = 1000, height = 800)
pielou_a
dev.off()

png(("aim2_age_bray.png"), width = 7.5, height = 5, units = "in", res = 300)
age.bray.ord.pp
dev.off()

png(("aim2_age_jaccard.png"), width = 7.5, height = 5, units = "in", res = 300)
age.jaccard.ord.pp
dev.off()

png(("aim2_age_wunifrac.png"), width = 7.5, height = 5, units = "in", res = 300)
age.wunifrac.ord.pp
dev.off()

png(("aim2_age_unifrac.png"), width = 7.5, height = 5, units = "in", res = 300)
age.unifrac.ord.pp
dev.off()

##### SUMMARY #####

## Alpha: kruskal-wallis and anova results
## Beta: permanova and post-hoc results

kruskal_sh
TukeyHSD(anova_sh_vs_age_log)

kruskal_ch
TukeyHSD(anova_ch_vs_age_log)

kruskal_pd
TukeyHSD(anova_pd_vs_age_log)

kruskal_pielou
TukeyHSD(anova_pielou_vs_age_log)

perm_bray
head(df_bray)

perm_jaccard
head(df_jaccard)

perm_wunifrac
head(df_wunifrac)

perm_unifrac
head(df_unifrac)
