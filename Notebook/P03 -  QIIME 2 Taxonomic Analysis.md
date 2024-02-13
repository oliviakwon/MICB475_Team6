# Taxonomic analysis - Align database 
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads team6_rep-seq.qza \
  --o-classification team6_dysautonomia_taxonomy.qza

qiime metadata tabulate \
  --m-input-file team6_dysautonomia_taxonomy.qza \
  --o-visualization team6_dysautonomia_taxonomy.qzv

# Taxonomy barplots
qiime taxa barplot \
  --i-table team6_dysautonomia_table.qza \
  --i-taxonomy team6_dysautonomia_taxonomy.qza \
  --m-metadata-file /mnt/datasets/project_2/dysautonomia/dysautonomia_metadata.tsv \
  --o-visualization team6_dysautonomia_taxa-bar-plots.qzv
