#!/bin/bash

# define input directory
input_dir="./fasta_v1"
csv_file="./results.csv"

# Clear the log file and initialize the CSV file with column names
> "$log_file"
echo "Filename,Precision,Recall,F1,Kmap_Consensus,Meme_Consensus,Overlap" > "$csv_file"

# iterate all fasta files
for fasta_file in "$input_dir"/*.fasta; do

    filename=$(basename "$fasta_file" .fasta)
    # mkdir
    meme_dir="./meme_v1/$filename"
    kmap_dir="./kmap_v1/$filename" 

    # Run the Python script, capture its output
    py_output=$(python cal_simlarity_2.py --kmap_out "$kmap_dir" --meme_out "$meme_dir")

    # Append filename to the Python script output
    csv_line="$filename,$py_output"

    echo "$filename,$py_output" >> "$csv_file"
done

echo "All processes completed."
