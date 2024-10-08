#------------- load libraries-------------
library(phyloseq)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tibble)
library(vegan)
library(PMCMRplus)
library(ggsignif)
library(ggpubr)
library(rstatix)
#install.packages("ggExtra")
#install.packages("devtools")
#devtools::install_github("daattali/ggExtra")
library("ggExtra")
library(RColorBrewer)
library(ggside)
library(aplot)


#------------- load data and subset-------------

load("FD_rare_only_separate.Rdata")

#duplicate phylosesq object
phyloseq_FD<- FD_rare_only_separate 

#extract the sample data
sample_phyloseq_FD <- sample_data(phyloseq_FD)

# Retain specific column
phyloseq_FD_subset <- sample_phyloseq_FD [, c("Sample.Name", "Age.New.Bin", "Cage.ID", "Experiment.Group", 
                                              "Genotype", "Mouse.ID", "Phenotype.score", 
                                              "FD.severity", "Sex", "Weight.grams")]

FD_phyloseq_subset <- FD_rare_only_separate

#------------- SET COLOR ------------
## Define colour palette

# view colourblind-friendly palettes
display.brewer.all(colorblindFriendly = TRUE)

# Set colour palette
colour <- brewer.pal(8, "Set2")  # Replace "Set2" with any desired palette
low_colour <- colour[1]
high_colour <- colour[2]

#> low_colour
#[1] "#66C2A5"
#> high_colour
#[1] "#FC8D62"

#------------- bray curtis - SEX -------------
#plot the order in Male, Female
sample_data(FD_phyloseq_subset)$Sex <- factor((sample_data(FD_phyloseq_subset)$Sex),
                                              levels = c("Male", "Female"))

#create ordination for bray curtis
sex.bray.ord = ordinate(FD_phyloseq_subset, "PCoA", "bray", na.rm = TRUE) 

#plot ordination
sex.bray.ord.p = plot_ordination(FD_phyloseq_subset,
                                 sex.bray.ord, 
                                 color = "Sex")
#title = "Bray-Curtis dissimilarity across mouse sexes") 

#format your plot
sex.bray.main = sex.bray.ord.p + #change this
  geom_point(size = 3)+ 
  scale_color_discrete(name = "Mouse sex")+
  scale_color_brewer(palette = "Set2")+
  theme_bw(base_size = 18,
           base_line_size =0) +
  stat_ellipse(type = "norm")+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
sex.bray.main


#preview marginal boxplot
ggMarginal(sex.bray.main, 
           type ="boxplot", groupColour = TRUE, groupFill = TRUE)

#see the individual cages in addition to the sex data
sex.cage.bray.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                 sex.bray.ord, #change this
                                 color = "Cage.ID",
                                 shape="Sex") #change this to the sex category
#title = "Bray-Curtis dissimilarity across mouse sexes") #change this
# Load 'ggside' package
library(ggside)
#format your plot
sex.cage.bray.main = sex.cage.bray.ord.p + #change this
  geom_point(size = 5)+ 
  theme_bw(base_size = 18,
           base_line_size =0) 
sex.cage.bray.main

# extract axis 1 and 2 values
bray.axis<- as.data.frame(sex.bray.ord$vectors[,1:2]) # note that the values would be same if you use the same metrics
colnames(bray.axis)<- gsub("Axis", "bray.axis", colnames(bray.axis))
#match sample with sample data
FD.subset = sample_data(FD_phyloseq_subset)
FD.subset.df <- data.frame(FD.subset)
#make rows into a column for both datasets
FD.subset.complete<- rownames_to_column(FD.subset.df,var = "raw-sample-id") #33 13
bray.axis.complete<- rownames_to_column(bray.axis,var = "raw-sample-id") #33 3
#then add axis 1 and 2 pcoa values derived from bray into the FD.subset.df
bray.fd.join<- dplyr::full_join(FD.subset.complete, bray.axis.complete, by = "raw-sample-id")


#perform wilcox rank sum test for pcoa axis 1 and 2
#axis 1 (p-value = 4.401e-05)
wilcox.bray.axis1<- wilcox.test(bray.fd.join$bray.axis.1 ~bray.fd.join$Sex)

wilcox.bray.axis1 <- bray.fd.join %>%
  wilcox_test(bray.axis.1 ~ Sex) %>%
  add_significance()
wilcox.bray.axis1

#axis 2 (p-value = 0.07408)
wilcox.bray.axis2 <- bray.fd.join %>%
  wilcox_test(bray.axis.2 ~ Sex) %>%
  add_significance()
wilcox.bray.axis2

#plot axis 1 values
wilcox.bray.axis1.1 <- wilcox.bray.axis1 %>% add_xy_position(x = "Sex")
wilcox.bray.axis1.1.p<- ggplot(data =bray.fd.join, aes(x = Sex, y = bray.axis.1, colour =Sex))+
  geom_boxplot()+
  #geom_point(size = 3)+
  scale_color_discrete()+
  scale_color_brewer(palette = "Set2")+
  geom_signif(comparisons = list (c("Male", "Female")),
              y_position = c(0.3),
              annotations = c("****"),
              color = "#000000",
              textsize =5) +
  #stat_pvalue_manual(
  # wilcox.bray.axis1.1, label = "{p.signif}",
  # vjust = +2, hjust =-1, bracket.nudge.x = 0.35) +
  theme_bw(base_size = 18,
           base_line_size =0)+
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  coord_flip() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
wilcox.bray.axis1.1.p


#plot axis 2 values
wilcox.bray.axis2.1 <- wilcox.bray.axis2 %>% add_xy_position(x = "Sex")
wilcox.bray.axis2.1.p<- ggplot(data =bray.fd.join, aes(x = Sex, y = bray.axis.2, colour =Sex))+
  geom_boxplot()+
  #geom_point(size = 3)+
  scale_color_discrete()+
  scale_color_brewer(palette = "Set2")+
  geom_signif(comparisons = list (c("Male", "Female")),
              y_position = c(0.3),
              annotations = c("ns"),
              color = "#000000",
              textsize =5) +
  #stat_pvalue_manual(
  # wilcox.bray.axis2.1, label = "{p.signif}",
  # vjust = -1, bracket.nudge.y = 0.1
  #) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme_bw(base_size = 18,
           base_line_size =0)+
  theme(
    legend.position = "none", 
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank())
wilcox.bray.axis2.1.p

# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
aim2.sex.bray.p <- sex.bray.main %>% 
  insert_top(wilcox.bray.axis1.1.p, height = 0.2) %>% 
  insert_right(wilcox.bray.axis2.1.p, width = 0.2)
aim2.sex.bray.p 

# perform permanova 
dist.bray <- vegdist(t(otu_table(FD_phyloseq_subset)), method="bray") # Bray-curtis

FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results
sex_perm_bray <- adonis2(dist.bray ~ Sex,data = FD_phyloseq_sam_df)
sex_perm_bray #0.001 ***, 0.14868

#------------- jaccard -sex-------------
#plot the order in Male, Female
sample_data(FD_phyloseq_subset)$Sex <- factor((sample_data(FD_phyloseq_subset)$Sex),
                                              levels = c("Male", "Female"))

#create ordination for bray curtis
sex.jaccard.ord = ordinate(FD_phyloseq_subset, "PCoA", "jaccard", na.rm = TRUE) 

#plot ordination
sex.jaccard.ord.p = plot_ordination(FD_phyloseq_subset, 
                                    sex.jaccard.ord, 
                                    color = "Sex") 
#title = "Bray-Curtis dissimilarity across mouse sexes") 
# Load 'ggside' package
library(ggside)
#format your plot
sex.jaccard.main = sex.jaccard.ord.p + #change this
  geom_point(size = 5)+ 
  scale_color_discrete(name = "Mouse sex")+
  scale_color_brewer(palette = "Set2")+
  theme_bw(base_size = 18,
           base_line_size =0) +
  stat_ellipse(type = "norm")+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
sex.jaccard.main


#preview marginal boxplot
ggMarginal(sex.jaccard.main, 
           type ="boxplot", groupColour = TRUE, groupFill = TRUE)


# extract axis 1 and 2 values
jaccard.axis<- as.data.frame(sex.jaccard.ord$vectors[,1:2]) # note that the values would be same if you use the same metrics
colnames(jaccard.axis)<- gsub("Axis", "jaccard.axis", colnames(jaccard.axis))

#match sample with sample data
FD.subset = sample_data(FD_phyloseq_subset)
FD.subset.df <- data.frame(FD.subset)

#make rows into a column for both datasets
FD.subset.complete<- rownames_to_column(FD.subset.df,var = "raw-sample-id") #33 13
jaccard.axis.complete<- rownames_to_column(jaccard.axis,var = "raw-sample-id") #33 3

#then add axis 1 and 2 pcoa values derived from bray into the FD.subset.df
jaccard.fd.join<- dplyr::full_join(FD.subset.complete, jaccard.axis.complete, by = "raw-sample-id")

# perform wilcox rank sum test
#axis 1 (p-value = 4.401e-05)
wilcox.jaccard.axis1<- wilcox.test(jaccard.fd.join$jaccard.axis.1 ~jaccard.fd.join$Sex)

wilcox.jaccard.axis1 <- jaccard.fd.join %>%
  wilcox_test(jaccard.axis.1 ~ Sex) %>%
  add_significance()
wilcox.jaccard.axis1

#axis 2 (p-value = 0.07408)
wilcox.jaccard.axis2 <- jaccard.fd.join %>%
  wilcox_test(jaccard.axis.2 ~ Sex) %>%
  add_significance()
wilcox.jaccard.axis2

#plot axis 1 values
wilcox.jaccard.axis1.1 <- wilcox.jaccard.axis1 %>% add_xy_position(x = "Sex")
wilcox.jaccard.axis1.1.p<- ggplot(data =jaccard.fd.join, aes(x = Sex, y = jaccard.axis.1, colour =Sex))+
  geom_boxplot()+
  #geom_point(size = 3)+
  scale_color_discrete()+
  scale_color_brewer(palette = "Set2")+
  geom_signif(comparisons = list (c("Male", "Female")),
              y_position = c(0.35),
              annotations = c("****"),
              color = "#000000",
              textsize =5) +
  #stat_pvalue_manual(
  # wilcox.bray.axis1.1, label = "{p.signif}",
  # vjust = +3, hjust =-2, bracket.nudge.x = 2) +
  theme_bw(base_size = 18,
           base_line_size =0)+
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  coord_flip() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
wilcox.jaccard.axis1.1.p

#plot axis 2 values
wilcox.jaccard.axis2.1 <- wilcox.jaccard.axis2 %>% add_xy_position(x = "Sex")
wilcox.jaccard.axis2.1.p<- ggplot(data =jaccard.fd.join, aes(x = Sex, y = jaccard.axis.2, colour =Sex))+
  geom_boxplot()+
  #geom_point(size = 3)+
  scale_color_discrete()+
  scale_color_brewer(palette = "Set2")+
  geom_signif(comparisons = list (c("Male", "Female")),
              y_position = c(0.35),
              annotations = c("ns"),
              color = "#000000",
              textsize =5) +
  #stat_pvalue_manual(
  # wilcox.bray.axis2.1, label = "{p.signif}",
  # vjust = -1, bracket.nudge.y = 0.1
  #) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme_bw(base_size = 18,
           base_line_size =0)+
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
wilcox.jaccard.axis2.1.p

# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
aim2.sex.jaccard.p <- sex.jaccard.main %>% 
  insert_top(wilcox.jaccard.axis1.1.p, height = 0.2) %>% 
  insert_right(wilcox.jaccard.axis2.1.p, width = 0.2)
aim2.sex.jaccard.p #done

# perform permanova 
dist.jaccard <- vegdist(t(otu_table(FD_phyloseq_subset)), method="jaccard") # Bray-curtis

FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results
sex_perm_jaccard <- adonis2(dist.jaccard ~ Sex,data = FD_phyloseq_sam_df)
sex_perm_jaccard #0.001 ***, 0.11


#------------- weighted unifrac - sex -------------
#plot the order in Male, Female
sample_data(FD_phyloseq_subset)$Sex <- factor((sample_data(FD_phyloseq_subset)$Sex),
                                              levels = c("Male", "Female"))

#create ordination for bray curtis
sex.wu.ord = ordinate(FD_phyloseq_subset, "PCoA", "unifrac", weighted =T) 

#plot ordination
sex.uw.ord.p = plot_ordination(FD_phyloseq_subset,
                               sex.wu.ord, 
                               color = "Sex") 
#title = "Bray-Curtis dissimilarity across mouse sexes") 
# Load 'ggside' package
library(ggside)
#format your plot
sex.uw.main = sex.uw.ord.p + #change this
  geom_point(size = 5)+ 
  scale_color_discrete(name = "Mouse sex")+
  scale_color_brewer(palette = "Set2")+
  theme_bw(base_size = 18,
           base_line_size =0) +
  stat_ellipse(type = "norm")+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
sex.uw.main

#preview marginal boxplot
ggMarginal(sex.uw.main, 
           type ="boxplot", groupColour = TRUE, groupFill = TRUE)


# extract axis 1 and 2 values
uw.axis<- as.data.frame(sex.wu.ord$vectors[,1:2]) # note that the values would be same if you use the same metrics
colnames(uw.axis)<- gsub("Axis", "uw.axis", colnames(uw.axis))
#match sample with sample data
FD.subset = sample_data(FD_phyloseq_subset)
FD.subset.df <- data.frame(FD.subset)
#make rows into a column for both datasets
FD.subset.complete<- rownames_to_column(FD.subset.df,var = "raw-sample-id") #33 13
uw.axis.complete<- rownames_to_column(uw.axis,var = "raw-sample-id") #33 3
#then add axis 1 and 2 pcoa values derived from bray into the FD.subset.df
uw.fd.join<- dplyr::full_join(FD.subset.complete, uw.axis.complete, by = "raw-sample-id")


#perform wilcox rank sum test
#axis 1 (p-value = 0.0000089 ****)
wilcox.uw.axis1<- wilcox.test(uw.fd.join$uw.axis.1 ~uw.fd.join$Sex)

wilcox.uw.axis1 <- uw.fd.join %>%
  wilcox_test(uw.axis.1 ~ Sex) %>%
  add_significance()
wilcox.uw.axis1

#axis 2 (p-value = 0.0819 ns)
wilcox.uw.axis2 <- uw.fd.join %>%
  wilcox_test(uw.axis.2 ~ Sex) %>%
  add_significance()
wilcox.uw.axis2

#plot axis 1 values

#wilcox.uw.axis1.1 <- wilcox.uw.axis1 %>% add_xy_position(x = "Sex")

wilcox.uw.axis1.1.p<- ggplot(data =uw.fd.join, aes(x = Sex, y = uw.axis.1, colour =Sex))+
  geom_boxplot()+
  #geom_point(size = 3)+
  #scale_color_discrete()+
  scale_color_brewer(palette = "Set2")+
  geom_signif(comparisons = list (c("Male", "Female")),
              y_position = c(0.1),
              annotations = c("***"),
              color = "#000000",
              textsize =5) +
  #stat_pvalue_manual(
  # wilcox.bray.axis1.1, label = "{p.signif}",
  # vjust = +5, hjust =-1, bracket.nudge.x = -1) +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme_bw(base_size = 18,
           base_line_size =0)+
  coord_flip() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
wilcox.uw.axis1.1.p 

#plot axis 2 values
wilcox.uw.axis2.1.p<- ggplot(data =uw.fd.join, aes(x = Sex, y = uw.axis.2, colour =Sex))+
  geom_boxplot()+
  #geom_point(size = 3)+
  scale_color_discrete()+
  scale_color_brewer(palette = "Set2")+
  geom_signif(comparisons = list (c("Male", "Female")),
              y_position = c(0.1),
              annotations = c("ns"),
              color = "#000000",
              textsize =5) +
  #stat_pvalue_manual(
  # wilcox.bray.axis1.1, label = "{p.signif}",
  # vjust = +5, hjust =-1, bracket.nudge.x = -1) +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme_bw(base_size = 18,
           base_line_size =0)+
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
wilcox.uw.axis2.1.p

# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
aim2.sex.uw.marginal.p <- sex.uw.main %>% 
  insert_top(wilcox.uw.axis1.1.p, height = 0.2) %>% 
  insert_right(wilcox.uw.axis2.1.p, width = 0.2)
aim2.sex.uw.marginal.p #done

# perform permanova 
dist.uw <- UniFrac(FD_phyloseq_subset, weighted=TRUE) # Weighted UniFrac

# view permanova results
perm_wunifrac <- adonis2(dist.uw ~ Sex, data=FD_phyloseq_sam_df)
perm_wunifrac #0.001 ***, 0.23



#------------ bray curtis- age --------------
#plot the order in young, middle, then old
sample_data(FD_phyloseq_subset)$Age.New.Bin <- factor((sample_data(FD_phyloseq_subset)$Age.New.Bin),
                                                      levels = c("young", "middle", "old"))

levels(sample_data(FD_phyloseq_subset)$Age.New.Bin) <- c("Young", "Middle", "Old")

#create ordination for bray curtis
age.bray.ord = ordinate(FD_phyloseq_subset, "PCoA", "bray", na.rm = TRUE) 

#plot ordination
age.bray.ord.p = plot_ordination(FD_phyloseq_subset, #change this
                                 age.bray.ord, #change this
                                 color = "Age.New.Bin") 
#title = "Bray-Curtis dissimilarity across mouse ages")

#format your plot
age.bray.main = age.bray.ord.p + #change this
  geom_point(size = 3)+ 
  # scale_color_discrete(name = "Mouse age")+
  scale_color_brewer(palette = "Set2", name = "Age")+
  theme_bw(base_size = 18,
           base_line_size =0) +
  stat_ellipse(type = "norm")+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
age.bray.main

# preview adding marginal plot 
ggMarginal(age.bray.main, 
           type ="boxplot", groupColour = TRUE, groupFill = TRUE)

# extract axis 1 and 2 values - same as sex for bray
bray.axis<- as.data.frame(age.bray.ord$vectors[,1:2]) # note that the values would be same if you use the same metrics
colnames(bray.axis)<- gsub("Axis", "bray.axis", colnames(bray.axis))
#match sample with sample data - same as sex for bray
FD.subset = sample_data(FD_phyloseq_subset)
FD.subset.df <- data.frame(FD.subset) 
#make rows into a column for both datasets - same as sex for bray
FD.subset.complete<- rownames_to_column(FD.subset.df,var = "raw-sample-id") #33 13
bray.axis.complete<- rownames_to_column(bray.axis,var = "raw-sample-id") #33 3
#then add axis 1 and 2 pcoa values derived from bray into the FD.subset.df - same as sex for bray
bray.fd.join<- dplyr::full_join(FD.subset.complete, bray.axis.complete, by = "raw-sample-id")

# perform permanova 
dist.bray <- vegdist(t(otu_table(FD_phyloseq_subset)), method="bray") # Bray-curtis

FD_phyloseq_sam = sample_data(FD_phyloseq_subset)
FD_phyloseq_sam_df <- data.frame(FD_phyloseq_sam)

# view permanova results (0.001 ***)

perm_bray <- adonis2(dist.bray ~ Age.New.Bin,data=FD_phyloseq_sam_df)

#plot axis 1 values -Kruskal-Wallis rank sum test (p-value = 0.01119 )
kruskal.bray.axis1<- kruskal.test(bray.axis.1 ~ Age.New.Bin, data = bray.fd.join)

#then use wilcox test
wilcox.bray.axis1 <- bray.fd.join %>%
  wilcox_test(bray.axis.1 ~ Age.New.Bin) %>%
  add_significance()
wilcox.bray.axis1

bray.age.axis1.p<- ggplot(data =bray.fd.join, aes(x = Age.New.Bin, y = bray.axis.1, colour =Age.New.Bin))+
  geom_boxplot()+
  #geom_point(size = 3)+
  #scale_color_discrete()+
  scale_color_brewer(palette = "Set2", , name = "Age")+
  geom_signif(comparisons = list(c("Old","Middle"), c("Young", "Middle"), c("Young","Old")),
              y_position = c(0.4, 0.3, 0.5),
              annotations = c("***","**","**"),
              color = "#000000") +
  #stat_pvalue_manual(
  # wilcox.bray.age.axis1, label = "{p.adj.signif}",
  #vjust = +2, hjust =-1, bracket.nudge.x = 0.35) +
  theme_bw(base_size = 18,
           base_line_size =0)+
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  coord_flip() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
bray.age.axis1.p

#plot axis 2 values - Kruskal-Wallis rank sum test (p-value = 0.1499), no need to do post hoc
kruskal.bray.axis2<- kruskal.test(bray.axis.2 ~ Age.New.Bin, data = bray.fd.join)

wilcox.bray.axis2 <- bray.fd.join %>%
  wilcox_test(bray.axis.2 ~ Age.New.Bin) %>%
  add_significance()
wilcox.bray.axis2

bray.age.axis2.p<- ggplot(data =bray.fd.join, aes(x = Age.New.Bin, y = bray.axis.2, colour =Age.New.Bin))+
  geom_boxplot()+
  #geom_point(size = 3)+
  #scale_color_discrete()+
  scale_color_brewer(palette = "Set2", name = "Age")+
  geom_signif(comparisons = list(c("Old","Middle"), c("Young", "Middle"), c("Young","Old")),
              y_position = c(0.4, 0.3, 0.5),
              annotations = c("ns","ns","ns"),
              color = "#000000") +
  #stat_pvalue_manual(
  # wilcox.bray.age.axis1, label = "{p.adj.signif}",
  #vjust = +2, hjust =-1, bracket.nudge.x = 0.35) +
  theme_bw(base_size = 18,
           base_line_size =0)+
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
bray.age.axis2.p

# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
library(aplot)
aim2.age.bray.marginal1.p <- age.bray.main %>% 
  insert_top(bray.age.axis1.p, height = 0.2) %>% 
  insert_right(bray.age.axis2.p, width = 0.2)
aim2.age.bray.marginal1.p
