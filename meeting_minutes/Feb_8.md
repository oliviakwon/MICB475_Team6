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

### Meeting Notes ###
* Complete human data is available → server
* Evelyn: 99 is not a bad sample size 
* Ask about the workflow we have rn (QIIME2 and R analysis) → is there anything we need to add/cut
* Retrain classifier?
    * Depends
    * Have to look at variable region
    * Can retrain or use the V4 that is used in the original paper
    * Avril: there is a pre-trained one on the server
    * It’s okay to use the parameters previous research one
    * Use “silva 815-330” from the server← check that the numbers are the same
* Refraction depth
    * Depend on aim 1 → dependent on variables (each variable would have different depth)
    * Could run refraction curve -> generate a rarefaction curve first to see if parameters or samples drop off
    * If the quality is good, we could keep the same across variables
        * Ujemi paper used 251 = sampling depth which is considered high quality 
    * Download raw reads from Qiime2 to R (unrarefied)
    * In Qiime2 just generate the taxonomy (will generate phylogenetic tree during this stage too), table, req-seq
    * If we do refraction differently, how do we compare? - purely for alpha and beta diversity
        * Don’t rarefy in qiime -> put in R loop with variables of interest 
    * Qiime: table, taxonomy (with phylogeny tree), rep-seqs
        * Could start demultiplexing right now 
    * Don’t need to rarefy for aim 3-4 -> will need raw reads as input
* Document scripts used for QIIME on github
* Have two folders for QIIME scripts and R scripts
    * .sh files → for bash scripting
    * Have a separate output folder for .qza, figures, etc
* Aim 1: Avril will provide us with more info too
* Aim 3: Differential abundance OR  indicator species
    * Avril: Could choose one because they’re similar analyses → keep both in the proposal
    * Indicator: list of predictive species 
* Aim 4: Functional analysis depends on variable we choose → keep in the proposal too
* If there are too many “NA” → go back to refseq → blast to get the most updated result
* Github: 
    * Can bold and italicize text, difficult to colour (latex(?) format)
    * Does not need to make it fancy

### Task over the coming week ###
* Find classifier
* Work with QIIME to start with R
* Generate refraction curve from QIIME (preferred) or R
