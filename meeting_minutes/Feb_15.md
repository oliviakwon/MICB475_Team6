### Main Questions for Feb 15th Meeting ###
* Which univariate analysis formula to use (identifying the x and y values)?
  - Currently have two x_values instead of one category -> Should we do a bivariate analysis?
* Resources related to univariate analysis or go over the basic steps

### Updates ###
* Finished generating refraction curve using QIIME2
  - Demultiplexing -> max read depth: 251bp, no trimming
  - Denoising/Clustering (DADA2) -> no change in total # samples 
  - Taxonomic analysis -> excluded mitochondria and chloroplast sequences
    - Used V4 classifier on server (/silva-138-99-515-806-nb-classifier.qza)
  - Generated two alpha rarefaction curves (with and without eukaryotic sequences)
    - With max depth 60,000
    
* Github Documentation Structure
  - Detailed description of each QIIME2 step documented under MICB475_Team6/Notebook 
  - Code lines, output qzv files, output figures under MICB475_Team6/QIIME2

### Next Steps ###
* Next step(a) - Aim 1
  - Find variables significantly correlated with FD by Univariate (or Bivariate?) Regression Analysis in R
  - Resources we found with examples:
    - [Univariate](https://medium.com/@nsaeedster/basic-data-analysis-in-r-part-i-univariate-analysis-52048d4283e8)
      - Boxplot/histogram/frequency chart
    - [Bivariate](https://www.statology.org/bivariate-analysis-in-r/)
      - Scatterplots/Correlation Coefficients, Simple linear regression

* Next step(b) - Write project Proposal 
  - Introduction & Background: FD, metadata category of our focus
  - Research Question: Which variable affects the microbiome composition in mice the most and how does it impact the microbial dynamics?
  - Experimental Aims & Rationale
  - Approach
    - Aim 1: Find significant variables correlating with FD
      - Univariate Regression analysis
    - Aim 2: Diversity analysis on the found variable(s)
      - Richness and evenness
    - Aim 3 (if time permits): Differential abundance analysis & Indicator species analysis
      - DESeq2, indicspecies package
    - Aim 4 (if time permits): Functional analysis [PICRUSt](https://picrust.github.io/picrust/)
      - Metabolic pathway comparison
* Tentative Timeline Draft
  - ![alt text](https://github.com/oliviakwon/MICB475_Team6/blob/main/meeting_minutes/micb_475_timeline.png)
