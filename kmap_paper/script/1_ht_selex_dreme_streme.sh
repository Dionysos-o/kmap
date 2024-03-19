#!/bin/bash

input_dir="./fasta_dem"

log_file="$input_dir/process_times_dreme_streme.log"

> "$log_file"


for fasta_file in "$input_dir"/*.fasta; do

    filename=$(basename "$fasta_file" .fasta)
    echo "processing: $fasta_file" | tee -a "$log_file"


    work_dir="$input_dir/$filename"
    mkdir -p "$work_dir"

    fasta-unique-names -r "$fasta_file"


    # 
    echo "Running dreme for $fasta_file..." | tee -a "$log_file"
    { time dreme -mink 5 -maxk 16 -m 1 -p "$fasta_file" -oc "$work_dir/dreme_out"; } 2>> "$log_file"
    { time ceqlogo -i1 "$work_dir/dreme_out/dreme.txt" -o "$work_dir/dreme_out/logo.png" -f PNG; } 2>> "$log_file"


    # 
    echo "Running streme for $fasta_file..." | tee -a "$log_file"
    { time streme -minw 5 -maxw 16 -nmotifs 1 -p "$fasta_file" -oc "$work_dir/streme_out"; } 2>> "$log_file"
    { time ceqlogo -i1 "$work_dir/streme_out/streme.html" -o "$work_dir/streme_out/logo.png" -f PNG; } 2>> "$log_file"




    done
done

echo "All processes are complete. Times and average p-values recorded in $log_file."
