Feb 23, 2024 

### Proposal
- Gene naming: Genes are usually lower-cased, but literature has been calling this gene “ELP1”. Should we use ELP1, Elp1, or elp1?
- Do we need to mention packages when calling different functions?

### Research objective
- Hypothesis for Aim 3 and Aim 4: How should we formulate our hypothesis? Rubric said to not write predictions.

### Experimental aims and rationale
Aim 2:
-bar graph (categorical variable(s)) and an alluvial plot (numerical variable(s))
-Confirm with relative abundance as alpha diversity 
-beta diversity - not using PERMANOVA but bolxplots to compare various levels within a variable
-Should I test the distribution using qqplot before conducting non-parametric tests?
-Confirm the meaning of each indices
-Passive vs active voice in this aim

### Approach: 
Aim 2, 2B and 2D:
-2B: In QIIME2, visualise the table-no-mitochondria-no-chloroplast.qza (Aim 1D) via QIIME2 visualization tool to visualize the interactive plot of each chosen predictors and identify the sampling depth for each predictor maximizing the retained reads and sample size.
-2D: In RStudio, rarefy each phyloseq subset based on their optimal sampling depths
-Perform Tukey post-hoc test for kruskal wallis test
-Spearman’s rank correlation with r-squared and p value together?
Aim 4 (Clarify inputs to PICRUSt2)

### Code running
- Diagnose code - why is code not working . 
