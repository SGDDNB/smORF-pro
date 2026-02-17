#!/bin/bash

log_file="fasta/protein/clustalo_skips.log"

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

count_sequences() {
  awk 'BEGIN{count=0} /^>/{count++} END{print count}' "$1"
}

count_nonempty_sequences() {
  awk '
    /^>/ {if (seqlen>0) good++; seqlen=0; next}
    {seqlen += length($0)}
    END {if (seqlen>0) good++; print good+0}
  ' "$1"
}

for file in fasta/protein/*.txt; do
  base_name=$(basename "$file" .txt)
  total_seqs=$(count_sequences "$file")
  nonempty_seqs=$(count_nonempty_sequences "$file")

  if [ "$total_seqs" -lt 2 ]; then
    echo "$(timestamp) SKIP ${base_name}: only ${total_seqs} sequence(s) found" >> "$log_file"
    continue
  fi

  if [ "$nonempty_seqs" -lt 2 ]; then
    echo "$(timestamp) SKIP ${base_name}: fewer than 2 non-empty sequences (${nonempty_seqs})" >> "$log_file"
    continue
  fi

  echo "$base_name"
  clustalo -i "$file" --force --distmat-out="fasta/protein/dist_mat/${base_name}.mat" --full -o "fasta/protein/msa/${base_name}.msa" --percent-id
done


