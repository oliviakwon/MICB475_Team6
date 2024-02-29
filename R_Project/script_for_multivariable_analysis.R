

# This script is designed to rank multiple variables in terms of their contribution to microbial diversity using a PERMANOVA and looking at the
# R- squared value and pvalue. The R-squared statistic will indicate how much any one variable contributes to driving diversity.

# I have left lines unfinished for you to complete! They appear like:

variable =   #

#Good luck :) Email me if you get stuck.

#set your working directory
#setwd("R_project")


#import libraries 
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)



#Load in phyloseq object after filtering and rarification:
#                                       Filter for Chloroplast/Mitochondrial sequences
#                                       Filter to remove ASV's with counts fewer than 5
#                                       Filter to remove samples with fewer than 100 reads
#                                       Rarify phyloseq object to chosen sampling depth
#Code for this can be found throughout the Canvas modules!
load("FD_rare.RData")


#Take out the ASV matrix from the phyloseq object and convert it to a dataframe.
ASV_table = otu_table(FD_rare)
ASV_table = as.data.frame(otu_table(FD_rare))#convert to data.frame 

#Take out the metadata from the phyloseq object and convert it to a dataframe.
meta = sample_data(FD_rare)
meta = as.data.frame(sample_data(FD_rare))#convert to data.frame 
meta_df <- data.frame(meta)

#ps = prune_taxa(taxa_names(meta)[1],meta)
#psmelt = psmelt(ps)
# not working -> gives error Error in otu_table(physeq) : 
#argument "taxa_are_rows" is missing, with no default


#Filter the metadata to remove columns that are not related to nutrients
column_names <- colnames(meta)
formatted_colnames <- paste0(colnames(meta), collapse = '", "')
cat(formatted_colnames)

# need to change into data.frame 

meta_filt <- meta_df [, c("Sample.Name", "Age.Months", "Age.Days", "Age.New.Bin", "Test.Age.New.Bin", 
                       "Age.Bin", "Cage.ID", "Experiment.Group", "Genotype", 
                       "Mouse.ID", "Mut.Control.Ratio", "Phenotype.score", 
                       "FD.severity", "Test.FD.severity", "Disease.Bin", "Sex", 
                       "Treatment.type", "Weight.grams")]

FD = colnames(meta_filt)

adonis.res = list()   #Build an empty list that will be filled up by the loop

#Create a loop to go over each variable in the metadata.
for (i in 1:length(FD)){ #COMPLETE the for loop argument. You need to to loop through as many variables that are present in "nutrients". Use a range, IE (1:?)
  print(i)#Printing i to keep track of loop progress
  meta_no_missing_data <- meta_filt[complete.cases(meta_filt[, i]), ]#Remove the rows in metadata that contain missing data for the i'th variable
  
  samples = rownames(meta_no_missing_data) #Create a vector of all the samples in the metadata after removing NA's
  ASV_mat = ASV_table[,samples] #Filter the ASV_table to remove the same individuals that were filtered out for NA values 
                                #This is important because we need the number of individuals represented in the ASV table and metadata to be the same.
  
  ASV_mat = t(ASV_mat +1) # Transpose the ASV matrix so that the sample IDs become rownames. This is required for  adonis2()
  dis = vegdist(ASV_mat, method = "bray",diag =T, upper = T) #Create the distance matrix based on a bray-curtis dissimilarity model
  adonis.res[[i]] = vegan::adonis2(as.formula(paste("dis~",FD[i],sep = "")), data = meta_no_missing_data) #This line runs the PERMANOVA test and puts it into the empty list we made above.
}


#Create an empty table to import the R-squared and pvalue.
result = matrix(NA, nrow = length(FD), ncol =2)

#Loop though each variable and generate a table that can be plotted.
for (i in 1:length(FD)){ #This for loop argument will be the same as on line 60
  
  result[i,1] = adonis.res[[i]][1,3] #Grab the R-squared statistic
  result[i,2] = adonis.res[[i]][1,5]#Grab the pvalue

}

rownames(result) = c(FD)   #Convert the rowmanes to variables 
colnames(result) = c("R2", "Pvalue")#Change the column names TO "R2" AND "Pvalue"
result = data.frame(result, stringsAsFactors = F) #Convert it to a data.frame (easiest to work with when plotting)
result$Padjust = p.adjust(result$Pvalue, method = "fdr") #Generate an adjusted pvalue to correct for the probability of false positives
result$factor =  #Create another column with variable names
View(result)
###############################PLOTTING

#Try and generate the plot yourself using ggplot.
#Here is a skeleton of the plotting code. I wrote "FILLOUT" where information needs to be added.

#Also, filter the results table to only include significant variables with a pvalue<0.05

result_filtered = subset(result, Padjust < 0.05)#Write solution here

#I use reorder() in the plot below. This is how you can look up what it does. 
?reorder() #However, the best way is to google it.

ggplot(data =result_filtered , aes(x = reorder(rownames(result_filtered), -R2),y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() + ylab("Adonis R2") + xlab("Variables")


#Chris 

# subset - single taxon -> melt
# filter first OTU 
# melt 
