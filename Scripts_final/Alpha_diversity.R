
#---------- load libraries----------
library(tidyverse)
library(phyloseq)
library(ggsignif)
library(picante)
library(ggh4x)
library("car")


#----------load data and subset----------

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

#---------- SHANNON ----------

#check distribution - normally distributed
ggplot(samp_dat_wdiv, aes(sample = Shannon)) + 
  geom_qq() +  # Creates a Q-Q plot
  geom_qq_line() +  # Adds a reference line to the Q-Q plot
  facet_grid(~ Sex) +  # Facet by 'Sex' to create separate plots for 'Female' and 'Male'
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +  # Label your axes
  theme_classic()

#run glm with sex and mouse age as covariables
glm.sh<- glm(Shannon ~ Sex+ Age.New.Bin, data =samp_dat_wdiv,
             family = gaussian())
#intercept(3.729384): the expected value when all predictors are at their reference levels  
#SexMale(0.087889, female is the reference level) : the expected difference in the Shannon index for males compared to females is 0.087889, holding age constant.
#AIC(18.87):  a measure of the relative quality of the statistical model for a given set of data. A lower AIC indicates a better model

#extract p value (p = 0.548, ns)
summary(glm.sh)

# graph shannon indexes to sex levels and annotate significance
shannon_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=Shannon)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  ylab('Shannon index') +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(4.4),
              annotations = c("ns"), textsize = 10) +
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(size =24),
        strip.text.x = element_text(size = 28),
        strip.text.y = element_text(size = 28),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28))
shannon_sex

#---------- CHAO1 #----------
#check distribution - normally distributed
library("car")
ggplot(samp_dat_wdiv, aes(sample = "Chao1")) + 
  geom_qq() +  # Creates a Q-Q plot
  geom_qq_line() +  # Adds a reference line to the Q-Q plot
  facet_grid(~ Sex) +  # Facet by 'Sex' to create separate plots for 'Female' and 'Male'
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +  # Label your axes
  theme_classic()

#run glm with sex and mouse age as covariables
glm.chao<- glm(Chao1 ~ Sex+ Age.New.Bin, data =samp_dat_wdiv,
               family = gaussian())
#intercept(135.824): the expected value when all predictors are at their reference levels  
#SexMale(9.482, female is the reference level): the expected difference in the Shannon index for males compared to females is 9.482, holding age constant.
#AIC(302.4):  a measure of the relative quality of the statistical model for a given set of data. A lower AIC indicates a better model

#extract p value (p = 0.380, ns)
summary(glm.chao)

# graph Chao1 indexes to sex levels and annotate significance
chao1_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=Chao1)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  ylab('Chao1 index') +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(180),
              annotations = c("ns"), textsize = 10) +
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(size =24),
        strip.text.x = element_text(size = 28),
        strip.text.y = element_text(size = 28),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28))
chao1_sex


#---------- FAITH PHYLOGENETIC DIVERSITY ----------
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

#run glm with sex and mouse age as covariables
glm.pd <- glm(PD ~ Sex + Age.New.Bin, data =samp_dat_wdiv,
             family = gaussian())
glm.pd
#intercept(12.6298): the expected value when all predictors are at their reference levels
#SexMale(0.087889, female is the reference level): the expected difference in the PD index for males compared to females is 9.482, holding age constant.
#AIC(107.9): a measure of the relative quality of the statistical model for a given set of data. A lower AIC indicates a better model

#extract p value (p = 0.134, ns)
summary(glm.pd)

# run Wilxocon test on pd indexes
#wilcox_pd <- wilcox.test(PD ~ Sex, data=samp_dat_wdiv, exact = FALSE)
#wilcox_pd

# graph pd indexes to sex levels and annotate significance
pd_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=PD)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  ylab("Faith's Phylogenetic Diversity index") +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(15.3),
              annotations = c("ns"), textsize = 10) +
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(size =24),
        strip.text.x = element_text(size = 28),
        strip.text.y = element_text(size = 28),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28)))
pd_sex


#---------- PIELOU'S EVENNESS ----------

# look for the total number of species in the subsetted data
total_species <- as.numeric(ntaxa(FD_rare_only_separate))

# calculate pielou's evenness index
samp_dat_wdiv$pielou <- samp_dat_wdiv$Shannon/log(total_species)

#run glm with sex and mouse age as covariables
glm.pielou <- glm(pielou ~ Sex + Age.New.Bin, data =samp_dat_wdiv,
              family = gaussian())
glm.pielou
#intercept(0.597631): the expected value when all predictors are at their reference levels
#SexMale(0.014084, female is the reference level): the expected difference in the PD index for males compared to females is 9.482, holding age constant.
#AIC(-102): a measure of the relative quality of the statistical model for a given set of data. A lower AIC indicates a better model

#extract p value (p = 0.534, ns)
summary(glm.pielou)

# run Wilxocon test on pielou's evenness indexes
#wilcox_pielou <- wilcox.test(pielou ~ Sex, data=samp_dat_wdiv, exact = FALSE)
#wilcox_pielou

# graph pd indexes to sex levels and annotate significance
pielou_sex <- ggplot(samp_dat_wdiv, aes(x=`Sex`, y=pielou)) +
  scale_x_discrete(limits =c("Male", "Female")) +
  geom_boxplot() +
  geom_point(size = 2) +
  ylab("Pielou's Evenness index") +
  xlab('Sex') +
  theme(axis.text.x=element_text(angle=, vjust=0.5, hjust=0)) +
  geom_signif(comparisons = list(c("Male","Female")),
              y_position = c(0.71),
              annotations = c("ns"), textsize = 10) +
 theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(size =24),
        strip.text.x = element_text(size = 28),
        strip.text.y = element_text(size = 28),
        axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28))
pielou_sex


#---------- EXPORT PLOTS AS PNGS ----------
png(("aim2_sex_shannon.glm.png"), width = 500, height = 500)
shannon_sex
dev.off()

png(("aim2_sex_chao1.glm.png"), width = 500, height = 500)
chao1_sex
dev.off()

png(("aim2_sex_pd.glm.png"), width = 500, height = 500)
pd_sex
dev.off()

png(("aim2_sex_pielou.glm.png"), width = 500, height = 500)
pielou_sex
dev.off()



#---------- SUMMARY: Wilcoxon results #----------
wilcox_sh

wilcox_ch

wilcox_pd

wilcox_pielou
