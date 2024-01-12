cat M632_S90_L001_R1_001.fastq | convertseq fastq fasta > M632_S90_L001_R1_001.fasta

seqstat M632_S90_L001_R1_001.fasta
keeplong -h
cat M632_S90_L001_R1_001.fasta | keeplong fasta 80 > M632_S90_L001_R1_001_80bp.fasta

cat M632_S90_L001_R1_001_80bp.fasta

seqstat  M632_S90_L001_R1_001_80bp.fasta > 



blastn -query M632_S90_L001_R1_001_80bp.fasta -db Database/Enterovirus_data -out M632_blastn_90id.txt -qcov_hsp_perc 95 -max_target_seqs 1 -outfmt 6 -num_threads 64 -perc_identity 90

