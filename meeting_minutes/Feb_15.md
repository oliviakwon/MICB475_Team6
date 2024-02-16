### Main Questions for Feb 15th Meeting ###
* Which univariate analysis formula to use (identifying the x and y values)?
  - Currently have two x_values instead of one category -> Should we do a bivariate analysis?
* Resources related to univariate analysis or go over the basic steps
* Questions regarding hw assignment 5 and 6

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
      - PERMANOVA
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
 
### Meeting notes ###
* Beta-diversity: in R
* Just use the no mitochondria no chloroplast data
* Parkenstein’s: univariate correlation is ok → keep it more simple
* Evelyn has script for univariate analysis → sent to us (this is our next step)
    * Similar to first R assignment
    * Make a list of variables → run the loop from Chris → come up with a table (p values) → plot (to support that there are no weird outliers) and significance, 4x3 beta-diversity, correct the q values
        * Need to correct for big comparisons
        * How will we correct for it? Does it correct itself? FDR correction (use mutate from tidyverse and then use padjust to perform FDR generating Q value)
        * Don't need bray-curtis 
    * Don’t need to interpret - bring to next meeting 
* With the significant p values, we can choose which variable to work with for aims 2 and 3
* Try running Chris’ script ourselves
    * Need to rarefy
    * Check sampling depth - email Evelyn with table and get the ‘ok’
    * There should be one rarefied object
    * Rationale of why we think a specific sampling depth
* Proposal: 
    * Rubric: check assignment page
    * Timeline: what we have is fine
    * For timeline: should also add manuscript draft time into it 
* Do not remove things from taxonomy
* Indicator species is randomized: use “set.seed()” in Rstudio before ISA to obtain the same results every time running the code
* *First meeting after reading break: Get aim 1 done*: Should probably start working on it during reading week 

