### Main questions for Feb 8 meeting ###
- Ask about the workflow we have right now (QIIME2 and R analysis) → is there anything we need to add/cut (check below)
- Do we need to re-train the phylogenetic data from the server? 
- Refraction: What rarefaction depth will you choose for your diversity analyses? How many samples are discarded at this rarefaction depth? Do we need it for every category?
- Can we color code on GitHub text files?


### Main research question ###
Which variable affects the microbiome composition in mice the most and how does it impact the microbial dynamics (regression analysis - univariant) 

### Tentative workflow ###
* "QIIME2 quality check":
   * Demultiplex 16s rRNA
   * Train classifier (SILVA) 
      * the mouse data is already in V3 and 4; compare the primers to see if the regions match up → if yes, we could use the phylogenetic data from server
      * V4: the 515F primer (GTGCCAGCMGCCGCGGTAA) and 806R primer (GGACTACHVGGGTWTCTAAT) were used for amplification 
   * We could reference everything from the UJEMI paper

* "R analysis": 
   * Aim 1: Univariate Regression analysis to find significant variables that correlate with FD
      * Y value = bray curtis → based on 16S
      * X value = any category 
      * P value → table (with all variables) and 1 correlation figure with variable the most significant 
   * Aim 2: Perform diversity analysis on the chosen variable(s)
      * Richness and evenness based on 
      * Diversity matrix
   * Aim 3 (if time permits): Differential abundance analysis (genus) & Indicator species analysis (species)
      * Method: using DESeq2 (for differential)
      * Method: indicspecies package (species analysis) 
   * Aim 4 (if time permits): Functional analysis
      * Method: PICRUSt: Phylogenetic Investigation of Communities by Reconstruction of Unobserved States 
         * https://picrust.github.io/picrust/
      * Metabolic pathway comparison
