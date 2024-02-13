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
7. If mitochondria or chloroplast taxons are present, remove these sequences and generate a new table. 
   
## Output files:
Path: /data/team6_dysautonomia_data
1. team6_dysautonomia_taxonomy.qza
2. team6_dysautonomia_taxonomy.qzv
3. team6_dysautonomia_taxa-bar-plots.qzv
4. team6_dysautonomia_table-no-mitochondria-no-chloroplast.qza

## Results: 

## Discussion:

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
