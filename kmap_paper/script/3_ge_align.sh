#!/bin/bash

REFERENCE_FILE="./ge_result/ref_seq.fasta"
mkdir -p alignments
for i in {1..6}
do
    CLUSTER_FILE="./ge_result/c${i}.fasta"
    mkdir -p "alignments/cluster_${i}"
    awk '/^>/{if(x>0) close(out); x++; out="alignments/cluster_'${i}'/seq_"x".fasta"}; {print > out}' $CLUSTER_FILE

    for SEQ_FILE in alignments/cluster_${i}/*.fasta; do

        BASE_NAME=$(basename $SEQ_FILE)
        muscle -in "$SEQ_FILE" -in2 "$REFERENCE_FILE" -out "alignments/cluster_${i}/${BASE_NAME%.fasta}_aligned.fasta"
    done
done

echo "Alignment completed."
