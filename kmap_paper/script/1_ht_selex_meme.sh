#!/bin/bash

# define input directory
input_dir="./fasta_v1"

# remove old log file
log_file="./process_times.log"
> "$log_file"

# get number of fasta files for progress bar
total_files=$(find "$input_dir" -name '*.fasta' | wc -l)
processed_files=0

# iterate all fasta files
for fasta_file in "$input_dir"/*.fasta; do
    # Progress bar
    processed_files=$((processed_files + 1))
    echo -ne "Processing progress: $processed_files / $total_files\r"

    # get name of sample 
    filename=$(basename "$fasta_file" .fasta)
    echo "processing: $fasta_file" | tee -a "$log_file"

    # mkdir
    work_dir="./meme_v1/$filename"
    mkdir -p "$work_dir"

    fasta-unique-names -r "$fasta_file"

    # run meme 
    echo "Starting meme for $fasta_file at $(date)" | tee -a "$log_file"
    start_time=$(date +%s)

    # Run meme and record time in a separate file in work_dir
    time_log="$work_dir/time_log.txt"
    { time meme -dna "$fasta_file" -minw 5 -maxw 16 -oc "$work_dir"; } 2> "$time_log"

    end_time=$(date +%s)
    echo "Finished meme for $fasta_file at $(date)" | tee -a "$log_file"
    echo "Total time for $fasta_file: $((end_time - start_time)) seconds" | tee -a "$log_file"

done

echo "All processes completed."
