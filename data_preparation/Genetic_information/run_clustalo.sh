#!/usr/bin/env bash
set -euo pipefail

# ======================================
# Usage:
#   ./run_clustalo.sh <BASE_DIR>
#
# Example:
#   ./run_clustalo.sh results/fasta/protein
# ======================================

# -------- input argument --------
BASE_DIR="${1:-}"
if [[ -z "${BASE_DIR}" ]]; then
    echo "ERROR: BASE_DIR not provided"
    echo "Usage: $0 <BASE_DIR>"
    exit 1
fi

# -------- derived paths --------
ALN_DIR="${BASE_DIR}/aln"
TMP_DIR="${TMPDIR:-/tmp}"
# --------------------------------

mkdir -p "${ALN_DIR}"

for f in "${BASE_DIR}"/*.txt; do
    base=$(basename "$f" .txt)

    # keep only sequences >=10 aa
    awk '
        /^>/ {
            if (seq && length(seq) >= 10) {
                print header
                print seq
            }
            header=$0
            seq=""
            next
        }
        { seq = seq $0 }
        END {
            if (seq && length(seq) >= 10) {
                print header
                print seq
            }
        }
    ' "$f" > "${TMP_DIR}/${base}.filtered.fa"

    n=$(grep -c "^>" "${TMP_DIR}/${base}.filtered.fa")

    # need at least 2 sequences
    if [ "$n" -lt 2 ]; then
        echo "SKIP (too few long seqs): $f"
        continue
    fi

    clustalo \
        -i "${TMP_DIR}/${base}.filtered.fa" \
        -o "${ALN_DIR}/${base}.aln" \
        --force \
        || echo "FAILED: $f"
done
