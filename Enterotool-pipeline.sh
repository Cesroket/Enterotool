#!/bin/bash

# Variables
stamp="$1"
Fastq_input1="$2"
Fastq_input2="$3"
parent_dir="./Enterotool-${stamp}" 
rawdata_dir="${parent_dir}/00.Rawdata"
stats_dir="${parent_dir}/01.stats"
results_dir="${parent_dir}/02.results"
interfasta1="${stamp}_cleaned_paired_1.fastq"
interfasta2="${stamp}_cleaned_paired_2.fastq"
merged_output="${results_dir}/${stamp}_merged.fastq"
stats_file_before="${stats_dir}/stats_raw.txt"
stats_file_after="${stats_dir}/stats_trimmed.txt"
Blast_output="${results_dir}/${stamp}_blastn_95id.txt"
Edition_output="${results_dir}/${stamp}_abundance_table.tsv"
output="${results_dir}/${stamp}_abundance.results.tsv"
png_output="${results_dir}/${stamp}_graph.png"

mkdir -p "$rawdata_dir" "$stats_dir" "$results_dir"

cp "$Fastq_input1" "$rawdata_dir"
cp "$Fastq_input2" "$rawdata_dir"

seqkit stats "$Fastq_input1" "$Fastq_input2" > "$stats_file_before"
echo "Statistics before trimming saved to $stats_file_before"

trimmomatic PE -threads 4 "$Fastq_input1" "$Fastq_input2" \
               ${stamp}_cleaned_paired_1.fastq ${stamp}_cleaned_unpaired_1.fastq \
               ${stamp}_cleaned_paired_2.fastq ${stamp}_cleaned_unpaired_2.fastq \
               ILLUMINACLIP:adapters.fa:2:30:10 \
               LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

seqkit stats "$interfasta1" "$interfasta2" > "$stats_file_after"
echo "Statistics after trimming saved to $stats_file_after"

vsearch --fastq_mergepairs "$interfasta1" \
        --reverse "$interfasta2" \
        --fastqout "$merged_output" --fastq_minovlen 10

if [ $? -ne 0 ]; then
    echo "Error during merging of reads." >&2
    exit 1
else
    echo "Reads successfully merged: $merged_output"
fi

vsearch --fastq_filter "$merged_output" --fastaout "${results_dir}/input.fasta"
vsearch --uchime_denovo "${results_dir}/input.fasta" \
        --nonchimeras "${results_dir}/output_nonchimeras.fasta" \
        --chimeras "${results_dir}/output_chimeras.fasta"

# BLAST 
blastn -query "${results_dir}/output_nonchimeras.fasta" -db Database/Enterovirus_data_clean \
       -out "$Blast_output" -qcov_hsp_perc 97 -max_target_seqs 1 -outfmt 6 -perc_identity 97

# parse blast output
awk '{counter[$2]++} END {for (line in counter) print line "\t" counter[line]}' "$Blast_output" | sed -e s'/|/\t/'g > "$Edition_output"

# Normalize via bash, rnanorm, was not available.
echo "Normalizing abundance table to CPM..."
total_count=$(awk 'NR > 1 {sum += $3} END {print sum}' "$Edition_output")

if [[ $total_count -eq 0 ]]; then
    echo "Total count is zero. Cannot normalize." >&2
    exit 1
fi

awk -v total="$total_count" 'BEGIN {FS=OFS="\t"} 
    NR == 1 {
        print $0, "CPM"
    } 
    NR > 1 {
        cpm = ($3 / total) * 1000000
        print $0, cpm
    }' "$Edition_output" > "$output"

echo "Normalization completed. Normalized output saved to $output"

#i think this is the best option, or at least, the essiest one .
python3 - <<EOF
import pandas as pd
import matplotlib.pyplot as plt

file_path = "$output"
output_file_name = "$png_output"
df = pd.read_csv(file_path, sep='\t')
categorias_relevantes = ['EnterovirusA', 'EnterovirusB', 'EnterovirusC', 'EnterovirusD', 'EnterovirusE']
df['Genus'] = df['Genus'].apply(lambda x: x if x in categorias_relevantes else 'Others')
grouped_data = df.groupby('Genus')['CPM'].sum()
normalized_data = grouped_data / grouped_data.sum()

ax = normalized_data.plot(kind='bar', figsize=(10, 6), color='tab:blue', legend=False)
plt.xlabel('Genus')
plt.ylabel('Proportion')
plt.title('Relative Proportions of Genus')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(output_file_name)
EOF

echo "Enterotool Finished"