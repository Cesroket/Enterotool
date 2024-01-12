#!/bin/bash
seqkit fq2fa  M632_S90_L001_R1_001.fastq | seqkit seq -m 80 > output.fasta

#Database
#makeblastdb -in Database_enterovirus_clean.fasta -dbtype nucl -out Enterovirus_data

#blast
blastn -query output.fasta -db Database/Enterovirus_data_clean -out M632_blastn_90id.txt -qcov_hsp_perc 95 -max_target_seqs 1 -outfmt 6  -perc_identity 90 #-num_threads 64

#Part 2: edition.
awk '{contador[$2]++} END {for (linea in contador) print linea "\t" contador[linea]}' "M632_blastn_90id.txt" | sed -e s'/|/\t/'g > M632_ent_abundance.tsv

cat headers.tsv M632_ent_abundance.tsv > M632_ent_abundance_headers.tsv