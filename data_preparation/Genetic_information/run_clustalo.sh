#!/bin/bash

# Iteration over text files

for file in fasta/protein/*.txt; do
   base_name=$(basename "$file" .txt)
   echo "$base_name"
   ./clustalo -i "$file" --distmat-out="fasta/protein/dist_mat/${base_name}.mat" --full -o "fasta/protein/msa/${base_name}.msa" --percent-id
done


