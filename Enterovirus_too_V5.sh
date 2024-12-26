#!/bin/bash

#Variables
stamp="$1"
Fastq_input1="$2"
Fastq_input2="$3"
interfasta1="${stamp}_cleaned_paired_1.fastq"
interfasta2="${stamp}_cleaned_paired_2.fastq"
output_dir="./${stamp}_output"
stats_file="${output_dir}/${stamp}_stats.txt"
Fasta_output="${output_dir}/spades.output/contigs.fasta"
Blast_output="${output_dir}/${stamp}_blastn_95id.txt"
Edition_output="${output_dir}/${stamp}_ent_abundance.tsv"
output="${output_dir}/${stamp}_abundance_output.txt"
png_output="${output_dir}/${stamp}_graph.png"

# Crear el directorio de salida si no existey
mkdir -p "$output_dir"

#trimmomatic
trimmomatic PE -threads 4 $Fastq_input1 $Fastq_input2 \
               ${stamp}_cleaned_paired_1.fastq ${stamp}_cleaned_unpaired_1.fastq \
               ${stamp}_cleaned_paired_2.fastq ${stamp}_cleaned_unpaired_2.fastq \
               ILLUMINACLIP:adapters.fa:2:30:10 \
               LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

#Sapdes, ni idea si utilizare spades al final o no. 
spades -1 ${stamp}_cleaned_paired_1.fastq -2 ${stamp}_cleaned_paired_2.fastq -o ./${output_dir}/spades.output

#Chimeras out
vsearch --fastq_filter input.fastq --fastaout input.fasta
vsearch --uchime_denovo input.fasta --nonchimeras output_nonchimeras.fasta --chimeras output_chimeras.fasta


#Paso 2, stats
seqkit stats "$Fasta_output" > "$stats_file"

echo "Fq2fa Done"

#Database
#makeblastdb -in Database_enterovirus_clean.fasta -dbtype nucl -out Enterovirus_data

#blast
blastn -query "$Fasta_output" -db Database/Enterovirus_data_clean -out "$Blast_output" -qcov_hsp_perc 97 -max_target_seqs 1 -outfmt 6  -perc_identity 97 

echo "Blast results save in "$Blast_output""

#Part 2: edition.
awk '{contador[$2]++} END {for (linea in contador) print linea "\t" contador[linea]}' "$Blast_output" | sed -e s'/|/\t/'g > "$Edition_output"

#normalization rnorm

echo -e "Genus\tSpecie\tCount\tCPM" > headers.tsv
cat headers.tsv "$Edition_output" > "$output"
rm headers.tsv

echo "Final results save in "$output""

rm *fastq

echo "Plotting ..."

python stack_graph.py "$output" "$png_output"

echo "Job Finish"