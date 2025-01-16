#!/bin/bash

# Variables
stamp="$1"
Fastq_input1="$2"
Fastq_input2="$3"
interfasta1="${stamp}_cleaned_paired_1.fastq"
interfasta2="${stamp}_cleaned_paired_2.fastq"
merged_output="${stamp}_merged.fastq" # Archivo fusionado
output_dir="./${stamp}_output"
stats_file="${output_dir}/${stamp}_stats.txt"
Fasta_output="${output_dir}/spades.output/contigs.fasta"
Blast_output="${output_dir}/${stamp}_blastn_95id.txt"
Edition_output="${output_dir}/${stamp}_ent_abundance.tsv"
output="${output_dir}/${stamp}_abundance_output.txt"
png_output="${output_dir}/${stamp}_graph.png"

# Crear el directorio de salida si no existe
mkdir -p "$output_dir"

# Trimmomatic: limpieza de lecturas
trimmomatic PE -threads 4 $Fastq_input1 $Fastq_input2 \
               ${stamp}_cleaned_paired_1.fastq ${stamp}_cleaned_unpaired_1.fastq \
               ${stamp}_cleaned_paired_2.fastq ${stamp}_cleaned_unpaired_2.fastq \
               ILLUMINACLIP:adapters.fa:2:30:10 \
               LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

# Mergear lecturas forward y reverse
vsearch --fastq_mergepairs ${stamp}_cleaned_paired_1.fastq \
        --reverse ${stamp}_cleaned_paired_2.fastq \
        --fastqout "$merged_output" --fastq_minovlen 10

# Verificar si el mergeado fue exitoso
if [ $? -ne 0 ]; then
    echo "Error during merging of reads." >&2
    exit 1
else
    echo "Reads merged successfully: $merged_output"
fi

# Chimeras out
vsearch --fastq_filter "$merged_output" --fastaout input.fasta
vsearch --uchime_denovo input.fasta --nonchimeras output_nonchimeras.fasta --chimeras output_chimeras.fasta

# Paso 2: Estadísticas
seqkit stats "$Fasta_output" > "$stats_file"

echo "Fq2fa Done"

# Base de datos BLAST
# makeblastdb -in Database_enterovirus_clean.fasta -dbtype nucl -out Enterovirus_data

# BLAST
blastn -query "$Fasta_output" -db Database/Enterovirus_data_clean -out "$Blast_output" -qcov_hsp_perc 97 -max_target_seqs 1 -outfmt 6  -perc_identity 97 

echo "Blast results saved in $Blast_output"

# Edición de resultados
awk '{contador[$2]++} END {for (linea in contador) print linea "\t" contador[linea]}' "$Blast_output" | sed -e s'/|/\t/'g > "$Edition_output"

# Normalización
echo -e "Genus\tSpecie\tCount\tCPM" > headers.tsv
cat headers.tsv "$Edition_output" > "$output"
rm headers.tsv

echo "Final results saved in $output"

# Limpiar archivos temporales
rm *fastq

echo "Plotting ..."

python stack_graph.py "$output" "$png_output"

echo "Job Finished"
