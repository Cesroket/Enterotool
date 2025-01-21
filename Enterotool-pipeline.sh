#!/bin/bash

# Variables
stamp="$1"
Fastq_input1="$2"
Fastq_input2="$3"
parent_dir="./Enterotool-${stamp}" 
rawdata_dir="${parent_dir}/00.Rawdata"
stats_dir="${parent_dir}/01.Stats"
results_dir="${parent_dir}/02.Results"
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
               LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:80

seqkit stats "$interfasta1" "$interfasta2" > "$stats_file_after"
echo "Statistics after trimming saved to $stats_file_after"

vsearch --fastq_mergepairs "$interfasta1" \
        --reverse "$interfasta2" \
        --fastqout "$merged_output" --fastq_minovlen 10

vsearch --fastq_filter "$merged_output" --fastaout "${results_dir}/input.fasta"
vsearch --uchime_denovo "${results_dir}/input.fasta" \
        --nonchimeras "${results_dir}/output_nonchimeras.fasta" \
        --chimeras "${results_dir}/output_chimeras.fasta"

# BLAST 
blastn -query "${results_dir}/output_nonchimeras.fasta" -db Database/Enterovirus_data_clean -out "$Blast_output" -qcov_hsp_perc 98 -max_target_seqs 1 -outfmt 6 -perc_identity 98

# parse blast output
awk '{counter[$2]++} END {for (line in counter) print line "\t" counter[line]}' "$Blast_output" | sed -e s'/|/\t/'g > "$Edition_output"


#headers-not working, i will fix this in v2
sed -i '1i Genogroups\tSerogroups*\tCount' "$Edition_output"

# Normalize via bash, rnanorm, was not available. 
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


#i think this is the best option, or at least, the essiest one.
python3 - <<EOF
import pandas as pd
import matplotlib.pyplot as plt
file_path = "$output"
output_file_name = "$png_output"
df = pd.read_csv(file_path, sep='\t')
cat_rel = ['EnterovirusA', 'EnterovirusB', 'EnterovirusC', 'EnterovirusD', 'EnterovirusE']
df['Genogroups'] = df['Genogroups'].apply(lambda x: x if x in cat_rel else 'Others')
grouped_data = df.groupby('Genogroups')['Count'].sum()
normalized_data = grouped_data / grouped_data.sum()
plot_data = pd.DataFrame({'Proportion': normalized_data})
ax = plot_data.T.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab10')
plt.legend(title='Genogroups', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('genogroups')
plt.ylabel('EV genogroups sequence abundances')
plt.title('$stamp')
plt.xticks([])
plt.tight_layout()
plt.savefig(output_file_name)
EOF

#cleaning, in v2 i will clear this.
mv ${stamp}_*_*paired_*.fastq  ${parent_dir}/00.Rawdata/

echo "Pipeline finished. Results saved to $results_dir"