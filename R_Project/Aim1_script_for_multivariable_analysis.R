# This script is designed to rank multiple variables in terms of their contribution to microbial diversity using a PERMANOVA and looking at the
# R- squared value and pvalue. The R-squared statistic will indicate how much any one variable contributes to driving diversity.


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
load("FD_rare_only_separate.RData")


#Take out the ASV matrix from the phyloseq object and convert it to a dataframe.
ASV_table = otu_table(FD_rare_only_separate)
ASV_table = as.data.frame(otu_table(FD_rare_only_separate))#convert to data.frame 

#Take out the metadata from the phyloseq object and convert it to a dataframe.
meta = sample_data(FD_rare_only_separate)
meta_df <- data.frame(meta)

#ps = prune_taxa(taxa_names(meta)[1],meta)
#psmelt = psmelt(ps)
# not working -> gives error Error in otu_table(physeq) : 
#argument "taxa_are_rows" is missing, with no default


#Filter the metadata to remove columns that are not related to FD

meta_filt <- meta_df [, c("Sample.Name", "Age.New.Bin", "Cage.ID", "Experiment.Group", 
                          "Genotype", "Mouse.ID", "Phenotype.score", 
                          "FD.severity", "Sex", "Weight.grams")]

FD = colnames(meta_filt)

adonis.res = list()   #Build an empty list that will be filled up by the loop

#Testing that Cage ID has an effect on beta diversity -> potential confounding variable. 
i = "Cage.ID"
meta_no_missing_data <- meta_filt[complete.cases(meta_filt[, i]), ]#Remove the rows in metadata that contain missing data for the i'th variable

samples = rownames(meta_no_missing_data) #Create a vector of all the samples in the metadata after removing NA's
ASV_mat = ASV_table[,samples] #Filter the ASV_table to remove the same individuals that were filtered out for NA values 
#This is important because we need the number of individuals represented in the ASV table and metadata to be the same.

ASV_mat = t(ASV_mat +1) # Transpose the ASV matrix so that the sample IDs become rownames. This is required for  adonis2()
dis = vegdist(ASV_mat, method = "bray",diag =T, upper = T) 
vegan::adonis2(as.formula(paste("dis~Cage.ID",sep = "")), data = meta_filt)

#Create a loop to go over each variable in the metadata.
for (i in 1:length(FD)){ #COMPLETE the for loop argument. You need to to loop through as many variables that are present in "nutrients". Use a range, IE (1:?)
  print(i)#Printing i to keep track of loop progress
  meta_no_missing_data <- meta_filt[complete.cases(meta_filt[, i]), ]#Remove the rows in metadata that contain missing data for the i'th variable
  
  samples = rownames(meta_no_missing_data) #Create a vector of all the samples in the metadata after removing NA's
  ASV_mat = ASV_table[,samples] #Filter the ASV_table to remove the same individuals that were filtered out for NA values 
                                #This is important because we need the number of individuals represented in the ASV table and metadata to be the same.
  
  ASV_mat = t(ASV_mat +1) # Transpose the ASV matrix so that the sample IDs become rownames. This is required for  adonis2()
  dis = vegdist(ASV_mat, method = "bray",diag =T, upper = T) #Create the distance matrix based on a bray-curtis dissimilarity model
  adonis.res[[i]] = vegan::adonis2(as.formula(paste("dis~",FD[i], "+Cage.ID", sep = "")), data = meta_no_missing_data) #This line runs the PERMANOVA test and puts it into the empty list we made above.
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
result$Factor =  c("Sample Name", "Age",
                   "Cage ID", "Experiment Group", "Genotype", "Mouse ID", "Phenotype Score",
                   "FD Severity", "Sex", "Weight (g)")#Create another column with variable names
View(result)

###############################PLOTTING

#Try and generate the plot yourself using ggplot.
#Here is a skeleton of the plotting code. I wrote "FILLOUT" where information needs to be added.

#Also, filter the results table to only include significant variables with a pvalue<0.05

result_filtered_Padjust = subset(result, Padjust < 0.05)#Write solution here
final_result <- result_filtered_Padjust[!(result_filtered_Padjust$Factor %in% c("Cage ID", "Experiment Group", "FD Severity")), ]

#I use reorder() in the plot below. This is how you can look up what it does. 
?reorder() #However, the best way is to google it.

result_bar_plot <- ggplot(data = final_result, aes(x = reorder(Factor, -R2), y = R2)) +
  geom_bar(stat = 'identity') +
  coord_flip() + 
  ylab("Adonis R2") + 
  xlab("Variables") +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 22), # Adjust title text size
    axis.text.x = element_text(size = 14),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),  # Adjust y-axis text size
    axis.title = element_text(size = 20),   # Adjust axis title text size
    plot.margin = margin(10, 10, 10, 10)    # Adjust bottom margin to make space for the title
  ) +
  ggtitle("Effect of Variables on Beta Diversity")

result_bar_plot


#####Saving######
ggsave(filename = "Aim1_result_bar_plot.png"
       , result_bar_plot
       , height=4, width=7)
save(result, file = "Aim1_result.Rdata")
save(result_filtered, file = "Aim1_result_filtered.Rdata")
