# P01 - QIIME 2 Demultiplexing: dysautonomia dataset 

Feb 10th, 2023

## Purpose:
To import and sort the 16S rRNA sequences based on barcode information in QIIME2 

## Material: 
1. QIIME2
2. dysautonomia_manifest.tsv (path: /mnt/datasets/project_2/dysautonomia)

## Method:
1. Connect to Cisco Anyconnect Secure VPN. Open the terminal and login to MICB 475 class server using the provided login credentials.
2. Create a dedicated directory for project 2: /data/team6_dysautonomia_data
3. Import and demultiplex the dysautonomia dataset using pair-end to create a qza file with
demultiplexed samples 
5. Create a visualization for the demultiplexed samples
6. Move the demux.qzv file to your local computer directory and view it using the QIIME 2

## Output files:
1. /data/team6_dysautonomia_data/team6_dysautonomia_paired_demux_seqs.qza 
2. /data/team6_dysautonomia_data/team6_dysautonomia_paired_demux_seqs.qzv - [file](/QIIME2/export/paired_demux_seqs.qzv) is also uploaded to this repository
   
## Results: 
The same result is generated from samples sequenced from both forward and reverse primers
1. Total number of reads: 6531021 
2. Total number of samples: 214
3. Range of sequencing depth: 7440-114131
4. Maximum read length (bp): 251
5. All the reads the same length of 251 bp
<img src="/QIIME2/figures/demultiplexed_sequence_counts_summary.png" height="500" width="500">
<img src="/QIIME2/figures/demultiplexed_sequence_length_summary.png" height="500" width="500">

## Discussion:
1. The maximum read depth (bp) was 251 while all 214 samples (the same number of samples using both forward and reverse primers) had 251 bp in length.
2. After demultiplexing, the medians of the quality score for all bases were consistently high with a similar median value of 30, presenting over 99.9% base call accuracy. Hence, this suggests that there would be no trimming required. The truncation length selected was 251 bp.
3. Samples sequenced using forward primers resulted in higher phred scores compared to that of the reverse primers, hence denoising step will use those amplified using forward primers

## Future direction:
1. Denoise sequences using the selected truncation length of 251 bp and determine ASVs with DADA2 (outputs: team6_dysautonomia_table.qza, denoising-stats.qza) 


