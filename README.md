# Enterotool

Enterotool is a fast bioinformatics pipeline for processing and analyzing sequencing data to identify and quantify Enterovirus genogroups.

This pipeline uses Seqkit for initial quality control (QC1) to check input data integrity. Trimmomatic then trims low-quality bases and short sequences. A second quality control step (QC2) validates the trimming. Vsearch merges forward and reverse reads, and removes chimeric sequences. The reads are aligned to a curated database using Blastn with strict thresholds (e.g., 98% nucleotide identity, default). An in-house bash script parses BLAST results, counts unique matches for each genogroup, and normalizes the data to Counts Per Million (CPM). Finally, Matplotlib visualizes the normalized data as a genogroups stacked bar chart, showing the relative abundances of identified genogroups.

This workflow was designed for rapid and efficient sequencing data processing, ensuring specific identification of EV genogroups. All analyses were performed on a Linux-based system (Ubuntu 22.04.5 LTS) to optimize computational performance and reproducibility.

![Sin título](https://github.com/user-attachments/assets/5094931f-be2f-4f83-aaf9-fca1707bd68f)

# Installing

To install all necessary dependencies and set up the environment, simply run the following command. It will guide you through the process:

```bash
chmod +x Installer.sh # Grant executable permissions
bash Installer.sh
```

# Usage

```bash
# Run the pipeline
bash Enterotool-pipeline.sh <Sample-code> <input-fastq_forward-R1.fastq> <input-fastq_reverse-R2.fastq>

```

- Ensure all dependencies are installed and available. If you are using Conda, simply run the provided installer for an easier setup.

- Verify that the database is located in the correct directory or update the database path in line 48 of the Enterotool-pipeline.sh script <-db Database/Enterovirus_data_clean>.

- You can use this pipeline with any other databases; for any other target, however, keep in mind that the specificity of the analysis depends on the specificity of the database used.

- If you encounter errors related to script execution (e.g., "command not found" or "bad interpreter"), run:

```bash
# Convert files to Unix format
dos2unix Enterotool-pipeline.sh
```

# Dependencies

The following tools and libraries are required to run the pipeline:

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
