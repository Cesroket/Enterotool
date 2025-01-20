# Enterotool

Enterotool is a fast bioinformatics pipeline for processing and analyzing sequencing data to identify and quantify Enterovirus genogroups.

This pipeline employs Seqkit for initial quality control (QC1) to ensure the integrity of the input data. Following this, the reads are processed using Trimmomatic to perform quality trimming, removing low-quality bases and sequences that do not meet the minimum read length or quality score thresholds. The trimmed reads undergo a second quality control step (QC2) to validate the effectiveness of the trimming process. The forward and reverse reads are merged using Vsearch, and potential chimeric sequences are identified and removed via a de novo approach. The processed reads are aligned against a curated database using Blastn, employing stringent thresholds (e.g., 98% nucleotide identity, by default) to ensure accurate genogroup identification. An in-house bash script parses the BLAST results, counting unique and valid matches for each genogroup. The resulting data is normalized to Counts Per Million (CPM) for sample consistency. Finally, the normalized data is visualized using Matplotlib, producing a genogroups stacked bar chart that highlights the relative abundances of identified genogroups.

This workflow was designed for rapid and efficient sequencing data processing, ensuring specific identification of EV genogroups. All analyses were performed on a Linux-based system (Ubuntu 22.04.5 LTS) to optimize computational performance and reproducibility.

![Sin título](https://github.com/user-attachments/assets/5094931f-be2f-4f83-aaf9-fca1707bd68f)

# Dependencies

The following tools and libraries are required to run the pipeline:
Bash.

Trimmomatic v0.39

Seqkit v2.9.0

Vsearch v2.29.1

Blastn v2.14.1

pandas v2.2.3

matplotlib v3.9.3

This database was obtained from the NCBI Virus repository (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), and manually curated.

# References

NCBI Virus repository: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/

Seqkit:
Wei Shen\*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.

Trimmomatic:
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics (Oxford, England), 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170

vsearch:
Rognes, T., Flouri, T., Nichols, B., Quince, C. & Mahé, F. (2016) 'VSEARCH: a versatile open source tool for metagenomics', PeerJ, 4, e2584. doi:10.7717/peerj.2584.

Blastn:
Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2.

Matplotlib:
J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.

# Licenses

This project is licensed under the GNU General Public License v3.0 (GPL-3.0).
