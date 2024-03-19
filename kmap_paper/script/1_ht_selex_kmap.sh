#!/bin/bash

# define input directory
input_dir="../../ht_selec_all/fasta_v1"

# remove old log file
log_file="../../ht_selec_all/process_times_kmap.log"
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
    work_dir="../../ht_selec_all/kmap_v1/$filename"
    mkdir -p "$work_dir"

    # run kmap 
    echo "Starting kmap for $fasta_file at $(date)" | tee -a "$log_file"
    start_time=$(date +%s)

    # Run kmap and record time in a separate file in work_dir
    time_log="$work_dir/time_log.txt"
    { time python kmap.py --function cons_id --input_fasta_file "$fasta_file" --min_k 5 --max_k 17 --visual_len 16 --out_dir "$work_dir"; } 2> "$time_log"


    end_time=$(date +%s)
    echo "Finished kmap for $fasta_file at $(date)" | tee -a "$log_file"
    echo "Total time for $fasta_file: $((end_time - start_time)) seconds" | tee -a "$log_file"

done

echo "All processes completed."
