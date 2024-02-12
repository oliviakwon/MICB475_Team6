# P01 - QIIME 2 Demultiplexing: dysautonomia dataset 

Feb 10th, 2023

## Purpose:
To import and sort the 16S rRNA sequences based on barcode information in QIIME2 

## Material: 
QIIME2 

## Method:
1. Connect to Cisco Anyconnect Secure VPN. Open the terminal and login to MICB 475 class server using the provided login credentials.
2. Create a dedicated directory for project 2: /data/team6_dysautonomia_data
3. Import and demultiplex the dysautonomia dataset using pair-end (path: /mnt/datasets/project_2/dysautonomia) to create a qza file with
demultiplexed samples (output path: data/team6_dysautonomia_data/team6_dysautonomia_paired_demux_seqs.qza).
5. Create a visualization for the demultiplexed samples (output: team6_dysautonomia_paired_demux_seqs.qzv).
6. Move the demux.qzv file to your local computer directory and view it using the QIIME 2

## Results:
### Demultiplexed sequence counts summary
<img src="Demultiplexed sequence counts summary.png" height="200" width="800" alt="Demultiplexed Sequence Counts Summary">

### Frequency Histogram

### Interactive quality plots

### Demultiplexed sequence length summary

## Discussion:
1. The maximum read depth (bp) was 251 while all 214 samples (the same number of samples using both forward and reverse primers) had 251 bp in length.
2. After demultiplexing, the medians of the quality score for all bases were consistently high with a similar median value of 30, presenting over 99.9% base call accuracy. Hence, this suggests that there would be no trimming required. The truncation length selected was 251 bp.




