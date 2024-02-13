# P03 - QIIME 2 Taxonomic Analysis and Taxnomoy-based filtering

Feb 12th, 2024

## Purpose:
To determine the taxonomic groups to which the ASVs correspond to, and to generate a taxonomy bar plot that illustrates the relative frequency of each taxonomic group. 
To perform taxonomy-based filtering and remove any mitochondria or chloroplast sequences. 

## Material: 
1. QIIME2
2. silva-138-99-515-806-nb-classifier.qza (/mnt/datasets/classifiers)
3. team6_req-seq.qza (/data/team6_dysautonomia_data)
4. team6_dysautonomia_table.qza (/data/team6_dysautonomia_data)
5. dysautonomia_metadata.tsv (/mnt/datasets/project_2/dysautonomia)

## Method:
1. Connect to Cisco Anyconnect Secure VPN. Open the terminal and login to MICB 475 class server using the provided login credentials.
2. Create a detached screen and name it "taxonomic analysis". 
3. Align the database.
4. Visualize the taxonomy by converting qza files to qzv.
5. Generate the taxonomy bar plot.
6. Transfer the taxonomy and taxonomy bar plot visualization files to local computer and view the taxon using view.QIIME2.org.
7. If mitochondria or chloroplast taxons are present, perform taxonomic-based filtering to remove these sequences and generate a new table.
8. Summarize the table without mitochondria and chloroplast.
9. Visualize the table without mitochondria and chloroplast by converting qza file to qzv.
10. Transfer the table without mitochondria and chloroplast to local computer and view the table using view.QIIME2.org. 
   
## Output files:
1. /data/team6_dysautonomia_datateam6_dysautonomia_taxonomy.qza
2. /data/team6_dysautonomia_datateam6_dysautonomia_taxonomy.qzv - [file](/QIIME2/export/taxonomy.qzv)
3. /data/team6_dysautonomia_datateam6_dysautonomia_taxa-bar-plots.qzv - [file](/QIIME2/export/taxa-bar-plots.qzv)
4. /data/team6_dysautonomia_datateam6_dysautonomia_table-no-mitochondria-no-chloroplast.qza
5. /data/team6_dysautonomia_datateam6_dysautonomia_table-no-mitochondria-no-chloroplast.qzv - [file](QIIME2/export/table-no-mitochondria-no-chloroplast.qzv)

## Results: 
* The taxonomy table had mitochondria taxon. 
* The resulting table without mitochondria and chloroplast sequences contains: 
   * Table summary: 
      * Number of samples: 214
      * Number of features: 2458
      * Total frequency: 3 123 176
   * Frequency per sample
      * Minimum frequency	3,164.0
      * 1st quartile	8,577.25
      * Median frequency	11,935.5
      * 3rd quartile	17,420.25
      * Maximum frequency	60,144.0
      * Mean frequency	14,594.280373831776
    
Taxonomic bar blot at level 7
> <img src="/QIIME2/figures/taxonomic_bar_plot_level_7.png" height="300">

Table summary and frequency per sample after filtering mitochondria and chloroplast
> <img src="/QIIME2/figures/table-no-mitochondria-no-chloroplast.png" height="500">

Attached images are generated from [https://view.qiime2.org/](https://view.qiime2.org/)

## Discussion:
As mitochondria taxon was found while examining the taxonomy table in view.QIIME2.org, taxonomic-based filtering was applied
to remove the mitochondria and chloroplast sequences and thus taxons from the table. 

## Future direction:
### Alpha Rarefaction 
#### Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences team6_rep-seq.qza \
  --o-alignment team6_dysautonomia_aligned-rep-seq.qza \
  --o-masked-alignment team6_dysautonomia_masked-aligned-rep-seq.qza \
  --o-tree team6_dysautonomia_unrooted-tree.qza \
  --o-rooted-tree team6_dysautonomia_rooted-tree.qza

#### Alpha-rarefaction: based on team6_dysautonomia_table-no-mitochondria-no-chloroplast.qza
qiime diversity alpha-rarefaction \
  --i-table  team6_dysautonomia_table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny team6_dysautonomia_rooted-tree.qza \
  --p-max-depth 60000 \
  --m-metadata-file /mnt/datasets/project_2/dysautonomia/dysautonomia_metadata.tsv \
  --o-visualization team6_dysautonomia_alpha-rarefaction-no-mitochondria-no-chloroplast.qzv

#### Alpha-rarefaction: based on team6_dysautonomia_table.qza
qiime diversity alpha-rarefaction \
  --i-table  team6_dysautonomia_table.qza \
  --i-phylogeny team6_dysautonomia_rooted-tree.qza \
  --p-max-depth 60000 \
  --m-metadata-file /mnt/datasets/project_2/dysautonomia/dysautonomia_metadata.tsv \
  --o-visualization team6_dysautonomia_alpha-rarefaction.qzv
