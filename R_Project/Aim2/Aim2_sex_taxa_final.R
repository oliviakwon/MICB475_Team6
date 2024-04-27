#---------- load libraries----------
library(tidyverse)
library(phyloseq)
library(ggsignif)
library(picante)
library(ggh4x)


#----------load data and subset----------

load("FD_final_only_separate.Rdata") #from aim 1

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
FD_otu_relative1<- rownames_to_column(FD_otu_relative, var ="ASV") #1105   60

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
ra.p<- ggplot(data =taxa_ra, aes(x= Sample.Name, y = Relative_Abundance))+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  ylab('Relative abundance (%)')+
  xlab("Sample name")+
  facet_nested(~ factor(Sex, , levels=c("Male", "Female")) + Genotype,scales = "free_x")+
  #facet_grid(~factor(Sex, , levels=c("Male", "Female")), scales = "free_x")+
  ggtitle("Relative abundance of phylum across mouse sexes")+
  theme_classic(base_size = 18)

ra.p


png(("taxa.plot.norarefy1.png"), width = 1500, height = 700)
ra.p
dev.off()
