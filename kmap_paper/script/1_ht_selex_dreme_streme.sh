#!/bin/bash

# 定义输入目录
input_dir="./fasta_dem"
# 定义日志文件
log_file="$input_dir/process_times_dreme_streme.log"

# 清空旧的日志文件
> "$log_file"

# 遍历文件夹中的所有.fasta文件
for fasta_file in "$input_dir"/*.fasta; do
    # 获取文件名，不包含路径
    filename=$(basename "$fasta_file" .fasta)
    echo "processing: $fasta_file" | tee -a "$log_file"

    # 创建一个以文件名命名的工作目录
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
