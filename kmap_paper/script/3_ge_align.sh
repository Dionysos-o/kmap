#!/bin/bash

# 参考序列文件
REFERENCE_FILE="./ge_result/ref_seq.fasta"

# 创建一个目录用于存储所有生成的文件
mkdir -p alignments

# 处理每个c1.fasta到c6.fasta文件
cat alignments/results_cluster_2/*_aligned.fasta > "alignments/cluster_2_all_aligned.fasta"
for i in {1..6}
do
    # 构建当前文件名
    CLUSTER_FILE="./ge_result/c${i}.fasta"
    # 创建一个目录来存储当前cluster的比对结果
    mkdir -p "alignments/cluster_${i}"

    # 用awk分离每个序列（假设每个序列以">"开头）
    awk '/^>/{if(x>0) close(out); x++; out="alignments/cluster_'${i}'/seq_"x".fasta"}; {print > out}' $CLUSTER_FILE

    # 对每个序列进行比对
    for SEQ_FILE in alignments/cluster_${i}/*.fasta; do
        # 获取文件基本名
        BASE_NAME=$(basename $SEQ_FILE)
        # 使用MUSCLE进行比对
        muscle -in "$SEQ_FILE" -in2 "$REFERENCE_FILE" -out "alignments/cluster_${i}/${BASE_NAME%.fasta}_aligned.fasta"
    done
done

echo "Alignment completed."
