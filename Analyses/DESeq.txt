
#### Loading packages #### 
library(DESeq2)
library(phyloseq) 
library(tidyverse)
library(RColorBrewer)

#### Load physoleq object ####
load('FD_final_only_separate.RData') 

# rename the phyloseq object
meta_FD <- FD_final_only_separate
# filter out control sample data
meta_FD_no_controls <- subset_samples(meta_FD, Genotype != "Control") 

# subset samples into Male vs Female mice
dysautonomia_sex <- subset_samples(meta_FD_no_controls, Sex %in% c("Male", "Female"))


#### DESeq for Male vs Female FD mice ####
# 1 was not added because no zero error 
deseq_dysautonomia_sex <- phyloseq_to_deseq2(dysautonomia_sex,~Sex)
DESEQ_sex <- DESeq(deseq_dysautonomia_sex)
res <- results(DESEQ_sex, tidy=TRUE, 
               #this will ensure that female is  reference group
               contrast = c("Sex","Male","Female"))
View(res)
# Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
sex_volcano_plot <- res %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(x = "Log2FoldChange", y = "-log10(padj")
ggsave(filename = "DESeq_volcanoplot_sex.png"
       , sex_volcano_plot
       , height=4, width=5) 

# Table of results
sigASVs_sex <- res %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_sex)
# Significant ASVs
sigASVs_sex_vec <- sigASVs_sex %>%
  pull(ASV)
# There is 10 ASV significantly different between the male and female groups



## Identifying the differentially abundant ASVs at the phylum level ##
sex_FD_DESeq <- prune_taxa(sigASVs_sex_vec,dysautonomia_sex)
sex_sigASVs <- tax_table(sex_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_sex) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
ggplot(sex_sigASVs) +
  geom_point(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle = 90))



## Identifying the differentially abundant ASVs at the genus level ##
# use this for export 
sex_FD_DESeq <- prune_taxa(sigASVs_sex_vec,dysautonomia_sex)
sex_sigASVs <- tax_table(sex_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_sex) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# Reformating the Genus and Phylum names to remove "p__"
sex_sigASVs$Phylum = gsub("p__","",sex_sigASVs$Phylum) 
sex_sigASVs$Genus = gsub("g__","",sex_sigASVs$Genus) 
sex_sigASVs_no_NA <- subset(sex_sigASVs, !is.na(Genus))
DESeq_sex <- ggplot(sex_sigASVs_no_NA) +
  scale_fill_brewer(palette = 'Set1') + 
  geom_bar(width = 0.5, size = 1, aes(x=reorder(Genus, log2FoldChange), 
                         y=log2FoldChange, fill= Phylum), stat="identity")+
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(size = 15, face = "bold"),
        axis.text.x.bottom = element_text(angle = -45)) +
  labs(x = "", y = "Log2FoldChange Male/Female") 

DESeq_sex
  
ggsave(filename = "DESeq_sex.png"
       , DESeq_sex
       , height=7, width=10) 



## Identifying the differentially abundant ASVs at the species level ##
sex_FD_DESeq <- prune_taxa(sigASVs_sex_vec,dysautonomia_sex)
sex_sigASVs <- tax_table(sex_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_sex) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))
ggplot(sex_sigASVs) +
  geom_point(aes(x=Species, y=log2FoldChange), stat="identity")+
  theme_bw()
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Log2FoldChange Female/Male")



