# P04 - QIIME 2 Alpha-rarefaction

Feb 12th, 2024

## Purpose:
To generate an alpha-rarefaction curve. 

## Material: 
1. QIIME2
2. team6_req-seq.qza (/data/team6_dysautonomia_data)
3. team6_dysautonomia_table-no-mitochondria-no-chloroplast.qza (/data/team6_dysautonomia_data)
4. team6_dysautonomia_table.qza (/data/team6_dysautonomia_data)
5. dysautonomia_metadata.tsv (/mnt/datasets/project_2/dysautonomia)

## Method:
1. Connect to Cisco Anyconnect Secure VPN. Open the terminal and login to MICB 475 class server using the provided login credentials.
2. Generate a tree for phylogenetic diversity analyses.
3. View the table to see maximum frequency in order to determine max depth for the alpha-rarefaction.
3. Generate alpha-rarefaction curve with a max depth of 60 000.
4. Transfer the alpha-rarefaction file to local computer and view the alpha-rarefaction curve using view.QIIME2.org.

## Output files:
1. /data/team6_dysautonomia_datateam6_dysautonomia_aligned-rep-seq.qza
2. /data/team6_dysautonomia_datateam6_dysautonomia_masked-aligned-rep-seq.qza
3. /data/team6_dysautonomia_datateam6_dysautonomia_unrooted-tree.qza
4. /data/team6_dysautonomia_datateam6_dysautonomia_rooted-tree.qza
5. /data/team6_dysautonomia_datateam6_dysautonomia_alpha-rarefaction.qzv - [file](/QIIME2/export/table-no-mitochondria-no-chloroplast.qzv)
6. /data/team6_dysautonomia_datateam6_dysautonomia_alpha-rarefaction-no-mitochondria-no-chloroplast.qzv - [file](/QIIME2/export/alpha-rarefaction-no-mitochondria-no-chloroplast.qzv)

## Results: 
Alpha-rarefaction graph
  <img src="/QIIME2/figures/alpha-rarefaction.png" height="400">

No mitochondria no chloroplast alpha-rarefaction graph 
  <img src="/QIIME2/figures/alpha-rarefaction-no-mitochondria-no-chloroplast.png" height="400">

Attached images are generated from [https://view.qiime2.org/](https://view.qiime2.org/)

## Discussion:
Two alpha rarefaction curves were generated, one using the table with taxonomic-based filtering in which mitochondria and
chloroplast sequences were filtered out and the other using the table without filtering applied. Both alpha-rarefaction
curves are similar. 

## Future direction:
Perform a univariate regression analysis to find significant variables that correlate with FD.
