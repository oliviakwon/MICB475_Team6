# Agenda

## 1. Aim 4: Revised graphs
* Control vs. FD - Heatmap (p < 0.05 and (Log2FoldChange > 1 or log2FoldChange < -1))
   * Ran DEseq2_function on metadata containing both genotypes
> <img src="/R_Project/Aim4/Control_FD_heatmap.png"> 

* FD - Heatmap (p < 0.05 and (Log2FoldChange > 1 or log2FoldChange < -1))
   * Ran DEseq2_function on metadata containing only mutant genotype
> <img src="/R_Project/Aim4/FD_heatmap.png"> 

* FD - Volcano plot (p < 0.05 and (Log2FoldChange > 1 or log2FoldChange < -1))
> <img src="/R_Project/Aim4/FD_volcano.png"> 
   
* FD - PCA plot (p < 0.05 and (Log2FoldChange > 1 or log2FoldChange < -1)):
> <img src="/R_Project/Aim4/FD_PCA.png"> 

* FD - log2FC Bar plot (p < 0.05 and (Log2FoldChange > 1 or log2FoldChange < -1))
> <img src="/R_Project/Aim4/FD_bar.png"> 

* FD - Error Bar (p < 0.05 and (Log2FoldChange > 1 or log2FoldChange < -1))
  * Error

> p <- pathway_errorbar(abundance = abundance_data_filtered %>% column_to_rownames("pathway"),
+                       daa_results_df = reduced_metacyc_daa_annotated_results_df,
+                       Group = metadata$Sex,
+                       ko_to_kegg = FALSE,
+                       p_values_threshold = 0.05,
+                       order = "group",
+                       select = NULL,
+                       p_value_bar = TRUE,
+                       colors = NULL,
+                       x_lab = "feature")
> p
Error in `guide_transform()`:
! <Guide> classes have been rewritten as <ggproto> classes.
The old S3 guide methods have been superseded.
Run `rlang::last_trace()` to see where the error occurred.


## 2. Aim 2 alpha diversity updated: check glm formula and interpretation
### Alpha diversity
* [Alpha diversity codes](/R_Project/Aim2/Aim2_sex_alpha_glm) with glm (please help to check if possible)
* Relative abundance plot with control and mutant separated
   > <img src="/R_Project/Aim2/aim2_sex_genotype_ra.png" height="300">
* Shannon
   > <img src="/R_Project/Aim2/aim2_sex_shannon.glm.png" height="300"> <img src="/R_Project/Aim2/aim2_sh_glm.stats.png" height="300">
* Chao1
   > <img src="/R_Project/Aim2/aim2_sex_chao1.glm.png" height="300"> <img src="/R_Project/Aim2/aim2_chao_glm.stats.png" height="300">
* Faith's PD
   > <img src="/R_Project/Aim2/aim2_sex_pd.glm.png" height="300"> <img src="/R_Project/Aim2/aim2_pd_glm.stats.png" height="300">
* Pielou's Evenness
   > <img src="/R_Project/Aim2/aim2_sex_pielou.glm.png" height="300"> <img src="/R_Project/Aim2/aim2_pielou_glm.stats.png" height="300">

### Beta diversity - check age stats with Avril
* Bray PCoA - Sex
   > <img src="/R_Project/Aim2/aim2.sex.bray.marginal.png" height="300"> <img src="/R_Project/Aim2/aim2_sex_bray_stats.png" height="300">
* Jaccard PCoA - Sex
  > <img src="/R_Project/Aim2/aim2.sex.jaccard.marginal.png" height="300"> <img src="/R_Project/Aim2/aim2_jaccard.stats.png" height="300">
* weighted unifrac PCoA - Sex
  > <img src="/R_Project/Aim2/aim2.sex.wunifrac.marginal.png" height="300"> <img src="/R_Project/Aim2/aim2_uw_stats.png" height="300">  
* Bray PCoA - Age
  > <img src="/R_Project/Aim2/aim2.age.bray.marginal1.png" height="300"> <img src="/R_Project/Aim2/aim2_age_bray_stats.png" height="300">

## 3. Figures on slides and selection for presentation
https://docs.google.com/presentation/d/1hx67dfx5R6RL3IdvGP3vc2pbf9by2Sjh4PAQTaD3kpI/edit#slide=id.g26b9945d9c0_0_5



