# Agenda for the day
## 1) Aim 4
* Dr. Sun confirmed that we should include the log2FC bar plot as a main figure 4B along with the PCA plot
* Relating taxonomy to metabolic pathways - couldn't find specific code example, and then creating a heatmap confusion
* Joined_df (contains ASV id, sample id, taxa info from tax_table(phyloseq_object)
> <img src="/R_Project/Aim4/Screen Shot 2024-04-11 at 12.27.14 PM.png"> 
* pathway abundance for each sample and pathways
> <img src="/R_Project/Aim4/Screen Shot 2024-04-11 at 12.27.28 PM.png"> 
* Problem: There are multiple taxa information for each sample. How can we match which taxa information is correct for the sample in the pathway abundance table? Or is there a different way to approach this?
* Question: Do we need to plot for just the 4 significant pathways in the log2FoldChange bar plot? Or for more pathways?

## 2) Manuscript
### Link to the manuscript
* https://docs.google.com/document/d/1byn9kZiTLz3-7zFFpsLjK64PreuiehcST-MeUw3RUJU/edit?usp=sharing
### Link to the supplementary 
* https://docs.google.com/document/d/1LUTa-eOtX0zQcXzSIzOPfXqJDsqNQYZHv3oeL_ZM4Sc/edit
### Questions
#### Introduction
* Unsure whether the evidence is strong enough to support the hypothesis (details are commented in the doc)
#### Methods
* Data processing using the QIIME2 pipeline, should I add that "a detailed script can be found in the github repository" in the methods section?
* Univariate regression analysis - should I be thanking Chris for providing the univariate regression analysis template? If so how should I be writing this?
* Not sure if alpha diversity should be included
#### Results
* Beta diversity - should it be a stand-alone part, or should it be mentioned with aim 1?
* Alpha diversity - not mentioned in main figures. Should we still include it as a stand-alone part in results?
* 1. Double check analysis for aim 1 age  (figure 1)
    1. whether I need to do post hoc with padjust
    2. The format in figure legend
#### Discussion
* Alpha diversity - unsure whether alpha diversity answers our research question
* Taxonomic composition:
    1. Formating the hypothesis from aim1 with previous studies
    2. Firmicutes vs Bacteroidota ratio: we have higher Firmicutes vs Bacteroidota ratio in female MICE: previously HUMAN studies show the opposite result between sex, and also higher F/B ratio = more prone to disease
#### Conclusion
* Is this an adequate conclusion? 
#### Limitations
* Go over all limitations we included - is there items that are missed?
* Is the highlighted limitation a future direction?
#### Future directions
* Go over all future directions we included - is there items that are missed?
#### Figures and supplementary
* Do we put all figures in a separate document, or do we include the main figures in the manuscript and put the supplementary figures into a separate one?
* figure 1. wording: covariates of sex and age
