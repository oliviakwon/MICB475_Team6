# Agenda: 

## 1. Aim 2
* Ran Bray-Curtis on all significant variables (i.e. age, sex, genotype and body weight)
* Based on the PCoA results, is it better we analyze age or sex?
   * However, there are 26 Females and 7 Males only
   * reminder: there are only 4 old mice if we use age category
* Based on the results of Aim 1 and 2, the rest of the aims be performed on the selected variable only right?

### Age
* PCoA
   > <img src="/R_Project/Aim2/aim2_age_bray_1.png" height="300">
* Permanova
   > <img src="/R_Project/Aim2/aim2_age_bray_perm.png" height="180">
### Sex
* PCoA
   > <img src="/R_Project/Aim2/aim2_sex_bray.png" height="300">
* Permanova
   > <img src="/R_Project/Aim2/aim2_sex_bray_perm.png" height="180">
### Genotype 
* PCoA
   > <img src="/R_Project/Aim2/aim2_genotype_bray.png" height="300">
* Permanova
   > <img src="/R_Project/Aim2/aim2_genotype_bray_perm.png" height="180">
### Weight 
* PCoA
   > <img src="/R_Project/Aim2/aim2_weight_bray.png" height="300">
* Permanova
   > <img src="/R_Project/Aim2/aim2_weight_bray_perm.png" height="180">

## 2. Aim 3
* Finished running ISA and DESeq, need to fix visualizing the volcano plots and bar plots.
* Do we need to present all data for the three categories of Age.New.Bin?
* make sure the subset and data manipulation was correct

## 3. Aim 4
* Finished running PICRUSt2
* Working on visulizing outputs with R
* Line of code not running, trying to figure it out (folliwing picrust_analysis.R on Canvas)
* 'sample_names = append(sample_names, "pathway")'


## 4. Figures for Aim 1 & Aim 2

What is the format requirement for figures? 

Aim 1: 
> <img src="/R_Project/Aim1_result_bar_plot.png" height="300">
Figure 1. Age is strongly associated with changes to the gut microbial composition in mouse models.  A univariate regression analysis was conducted on the mouse metadata. Age (p = 0.002) ranks highest among biological factors influencing gut microbial composition, followed by sex (p = 0.002), weight (g) (p = 0.003), and genotype (p = 0.004), in descending order of significance. 

Aim 2: 
#### Alpha diversity - Age
* Shannon
  > <img src="/R_Project/Aim2/aim2_age_shannon.png" height="300">
* Chao1
  > <img src="/R_Project/Aim2/aim2_age_chao1.png" height="300">
* Pielou's evenness
  > <img src="/R_Project/Aim2/aim2_age_pielou.png" height="300">
* Faith's phylogenetic diversity
  > <img src="/R_Project/Aim2/aim2_age_pd.png" height="300">

#### Beta diversity - Age
* Bray
  > <img src="/R_Project/Aim2/aim2_age_bray_1.png" height="300">
* Jaccard
  > <img src="/R_Project/Aim2/aim2_age_jaccard.png" height="300">
* Unweighted unifrac
  > <img src="/R_Project/Aim2/aim2_age_unifrac.png" height="300">
* Weighted unifrac
  > <img src="/R_Project/Aim2/aim2_age_wunifrac.png" height="300">

