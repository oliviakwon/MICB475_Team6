# P02 - QIIME 2 Denoising and Clustering Sequences into ASVs

Feb 12th, 2023

## Purpose:
To detect and correct sequencing errors, and to group the sequences into respective ASVs

## Material: 
1. QIIME2 
2. team6_dysautonomia_paired_demux_seqs.qza (/data/team6_dysautonomia_data)
3. dysautonomia_metadata.tsv (/mnt/datasets/project_2/dysautonomia)

## Method:
1. Connect to Cisco Anyconnect Secure VPN. Open the terminal and login to MICB 475 class server using the provided login credentials.
2. Create a detached screen and name it "denoising". 
3. Denoise and cluster the demultiplexed sequences using DADA2
4. Visualize the ASVs by converting qza files to qzv.
5. Transfer the visualization files to local computer and view the representative sequences and table.qzv using view.QIIME2.org

## Output files:
Path: /data/team6_dysautonomia_data
1. team6_dysautonomia_table.qza:
2. team6_dysautonomia_table.qzv
3. team6_rep-seqs.qza
4. team6_rep-seqs.qzv
5. denoising-stats.qza

## Results: AFTER DENOISING/CLUSTERING:
1. Total number of reads retained: 3,134,413  
2. Total number of ASVs: 2,463
3. Total number of samples: 214
4. Range of sequencing depth: 3,164 -60,202

## Discussion:
The total number of samples did not change since all the samples had reads of 150 bp in length and every read with 150 bp was retained.


## Future direction:
The clustered ASVs can be filtered to exclude any mitochondrial or chloroplast sequences, which will remove any eukaryotic data. 

