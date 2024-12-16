# Enterotool
Bioinformatic analyses were performed to identify the EV genogroups. A comprehensive database was constructed that included the full genomic variability of the target group obtained from the NCBI Virus repository (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). Quality statistics were obtained using Seqkit (Wei Shen et al., 2024). Raw sequencing reads were processed and quality-trimmed using Trimmomatic (v0.39), applying a minimum read length threshold of 80 base pairs and ensuring a minimum Phred quality score of Q30 (Bolger et al., 2014). Chimera assessment was performed using vsearch (2.29.1), excluding potential chimeric sequences via --uchime_denovo (Rognes et al., 2016). Genogroup identification and annotation were performed using the Blastn (Altschul et al., 1990). An in-house bash script, EnteroTool, was developed for automated quantification of genogroup abundance, parsing Blastn output files to count matches identified per sample. Unique matches, where a read aligns to a single genogroup without ambiguity, were considered in the quantification. Data normalization to counts per million (CPM) was performed using the rnanorm package (v2.1.0) (Zmrzlikar et al., 2023).  The workflow was designed to provide rapid and efficient processing of sequencing data. All analyses were performed on a Linux-based system (Ubuntu 22.04.5 LTS) to optimize computational performance and reproducibility.


#Bibliography
NCBI Virus repository: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/
Seqkit: Wei Shen et al., 2024.
Trimmomatic: Bolger et al., 2014.
vsearch: Rognes et al., 2016.
Blastn: Altschul et al., 1990.
RNAnorm package: Zmrzlikar et al., 2023.
Matplotlib: ...
