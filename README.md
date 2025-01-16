# Enterotool

This database was obtained from the NCBI Virus repository (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), and manually curated.

Raw sequencing reads were processed and quality-trimmed using Trimmomatic (v0.39), applying a minimum read length threshold of 80 base pairs and ensuring a minimum Phred quality score of Q30 (Bolger et al., 2014). Chimera assessment was performed using vsearch (2.29.1), excluding potential chimeric sequences via --uchime_denovo (Rognes et al., 2016). Genogroup identification and annotation were dependent on the 5'UTR region. Blastn (Altschul et al., 1990) was employed to align reads against the 5'UTR-based curated database with a threshold of 98% nucleotide sequence identity to ensure high-confidence matches. An in-house bash script was developed for automated quantification of genogroup abundance, parsing Blastn output files to count valid matches identified per sample. Unique matches, where a read aligned to a single genogroup without ambiguity and met the identity criterion, were included in the analysis. Data normalization to counts per million (CPM) was performed using bash.

The workflow was designed for rapid and efficient processing of sequencing data, ensuring specific identification of EV genogroups. All analyses were performed on a Linux-based system (Ubuntu 22.04.5 LTS) to optimize computational performance and reproducibility.

![Sin título](https://github.com/user-attachments/assets/5094931f-be2f-4f83-aaf9-fca1707bd68f)


# References

NCBI Virus repository: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/
Seqkit: Wei Shen et al., 2024.
Trimmomatic: Bolger et al., 2014.
vsearch: Rognes et al., 2016.
Blastn: Altschul et al., 1990.
Matplotlib: ...
