#!/usr/bin/env python3
"""
bio_parse_orf_gtf_no_length_restriction.py
=========================================

Purpose
-------
This script performs evolutionary conservation analysis of smORFs / iORFs
using multi-species whole-genome alignments (MAF files).

Given:
- a smORF-pro formatted GTF file (iORF / orfCDS features),
- UCSC multi-species MAF files (e.g. 100-way vertebrate alignment),

the script:
1. Extracts spliced nucleotide alignments for each iORF across species
2. Writes nucleotide FASTA and MAF alignments
3. Translates in-frame CDS sequences into protein FASTA files (one per iORF)
4. Computes basic conservation metrics:
   - start codon conservation
   - stop codon conservation
   - frame preservation
   - nucleotide percent identity vs reference genome

The resulting protein FASTA files are intended to be used as input for
Clustal Omega (see run_clustalo.sh).

Important notes
---------------
- No length restriction is applied to ORFs (intended for smORFs).
- MAF files are indexed on first run; this may take hours to days.
- Only in-frame orthologs are translated into proteins.
- This script does NOT run multiple sequence alignment itself.

Expected GTF format
-------------------
Feature types:
- iORF
- orfCDS

Required attributes:
- ORF_id
- iORF_id

Output structure
----------------
<outdir>/
  ├── fasta/
  │   ├── <iORF>.fa          (nucleotide alignments)
  │   └── protein/
  │       └── <iORF>.txt     (protein FASTA, multi-species)
  ├── maf/
  │   └── <iORF>.maf
  ├── start_codon_conservation.txt
  ├── stop_codon_conservation.txt
  ├── inframe.txt
  └── percentage_ids.txt

Example usage
-------------
python bio_parse_orf_gtf_no_length_restriction.py \
  --gtf Collapsed_hubner_fixed_V2.gtf \
  --maf-dir /data/ucsc/hg38/multiz100way/maf \
  --species species_names.txt \
  --outdir Conservation \
  --ref-genome hg38

"""

import os
import re
import argparse
import functools
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq


# --------------------------------------------------
# CLI
# --------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract smORF conservation from multi-species MAF alignments"
    )
    parser.add_argument("--gtf", required=True, help="smORF GTF file")
    parser.add_argument("--maf-dir", required=True, help="Directory with per-chromosome MAF files")
    parser.add_argument("--species", required=True, help="Species names table (GENOME / COMMON)")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current directory)")
    parser.add_argument("--ref-genome", default="hg38", help="Reference genome name in MAF (default: hg38)")
    return parser.parse_args()


# --------------------------------------------------
# GTF parsing
# --------------------------------------------------

def parse_gtf(fn):
    orfs = {}
    iorfs = {}
    orf_intervals = {}

    with open(fn) as gtf:
        for l in gtf:
            fields = l.rstrip().split("\t")
            interval = (fields[0], str(int(fields[3]) - 1), fields[4], fields[6])

            # FIXED attribute parsing
            attributes_split = {}
            for attr in fields[8].strip().rstrip(";").split("; "):
                key, value = attr.split(" ", 1)
                attributes_split[key] = value.strip('"')

            if "ORF_id" not in attributes_split:
                raise ValueError(f"Missing ORF_id in GTF line:\n{l}")

            if fields[2] == "iORF":
                orfs[attributes_split["ORF_id"]] = set()
                orf_intervals[attributes_split["ORF_id"]] = interval

            elif fields[2] == "orfCDS":
                if "iORF_id" not in attributes_split:
                    raise ValueError(f"Missing iORF_id in GTF line:\n{l}")

                orfs[attributes_split["ORF_id"]].add(attributes_split["iORF_id"])
                iorfs.setdefault(attributes_split["iORF_id"], []).append(interval)

    for iorf_id in iorfs:
        iorfs[iorf_id].sort(key=functools.cmp_to_key(sort_intervals))

    print(len(orfs))
    print(len(iorfs))
    return orfs, iorfs, orf_intervals



def sort_intervals(x, y):
    if x[3] == "+":
        return (x[1] > y[1]) - (x[1] < y[1])
    else:
        return (x[2] > y[2]) - (x[2] < y[2])


# --------------------------------------------------
# MAF extraction
# --------------------------------------------------

def get_maf(iorfs, maf_dir, ref_genome, fasta_dir, maf_outdir):
    all_sequences = {}
    idxes = {}

    chromosomes = [
        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chr20","chr21","chr22","chrX","chrY","chrM"
    ]

    # --- Use absolute paths for mafindex ---
    for chr in chromosomes:
        maf_path = os.path.join(maf_dir, f"{chr}.maf")
        idx_path = os.path.join(maf_dir, f"{chr}.mafindex")

        if not os.path.exists(idx_path):
            raise FileNotFoundError(f"Missing MAF index: {idx_path}")

        idxes[chr] = AlignIO.MafIO.MafIndex(
            idx_path,
            maf_path,
            f"{ref_genome}.{chr}"
        )

    # --- Extract spliced alignments ---
    for iorf_id, cds_intervals in iorfs.items():
        chrom = cds_intervals[0][0]   # e.g. "1", "X"
        idx = idxes[f"chr{chrom}"]
        strand = 1 if cds_intervals[0][3] == "+" else -1

        splices = idx.get_spliced(
            [int(i[1]) for i in cds_intervals],
            [int(i[2]) for i in cds_intervals],
            strand=strand
        )

        all_sequences[iorf_id] = {
            ot.id.split(".")[0]: ot.seq for ot in splices
        }

        AlignIO.write(splices,
                      os.path.join(fasta_dir, f"{iorf_id}.fa"),
                      "fasta")

        AlignIO.write(splices,
                      os.path.join(maf_outdir, f"{iorf_id}.maf"),
                      "maf")

    return all_sequences



# --------------------------------------------------
# Conservation metrics
# --------------------------------------------------

def process_start_codon(sequences, species):
    df = pd.DataFrame(index=sequences.keys(), columns=species["GENOME"])
    for iorf_id, seqs in sequences.items():
        if "hg38" not in seqs:
            continue
        human = str(seqs["hg38"]).lower()
        idx = [m.start() for m in re.finditer("[atgcn]", human)][:3]
        for sp, seq in seqs.items():
            df.at[iorf_id, sp] = "".join(seq.lower()[i] for i in idx)
    return df


def process_stop_codon(sequences, species):
    df = pd.DataFrame(index=sequences.keys(), columns=species["GENOME"])
    for iorf_id, seqs in sequences.items():
        if "hg38" not in seqs:
            continue
        human = str(seqs["hg38"]).lower()
        idx = [m.start() for m in re.finditer("[atgcn]", human)][-3:]
        for sp, seq in seqs.items():
            df.at[iorf_id, sp] = "".join(seq.lower()[i] for i in idx)
    return df


def process_sequences_inframe(sequences, species):
    df = pd.DataFrame(index=sequences.keys(), columns=species["GENOME"])
    for iorf_id, seqs in sequences.items():
        for sp, seq in seqs.items():
            df.at[iorf_id, sp] = len(seq.replace("-", "")) % 3 == 0
    return df


def process_sequences_pid(sequences, species, target="hg38"):
    df = pd.DataFrame(index=sequences.keys(), columns=species["GENOME"])
    for iorf_id, seqs in sequences.items():
        if target not in seqs:
            continue
        ref = seqs[target].lower()
        for sp, seq in seqs.items():
            if len(seq) != len(ref):
                continue
            match = total = 0
            for a, b in zip(ref, seq.lower()):
                if a == "-" or b == "-":
                    continue
                total += 1
                if a == b:
                    match += 1
            if total > 0:
                df.at[iorf_id, sp] = 100 * match / total
    return df


def export_protein_fastas(sequences, protein_dir):
    for iorf_id, seqs in sequences.items():
        out_path = os.path.join(protein_dir, f"{iorf_id}.txt")

        records = []
        for sp, seq in seqs.items():
            ungapped = seq.replace("-", "")
            if len(ungapped) % 3 != 0:
                continue

            protein = ungapped.translate(to_stop=True)

            # critical filter
            if len(protein) == 0:
                continue

            records.append((sp, protein))

        # write only if >= 2 sequences (Clustal Omega requirement)
        if len(records) < 2:
            continue

        with open(out_path, "w") as fh:
            for sp, protein in records:
                fh.write(f">{sp}\n{protein}\n")



# --------------------------------------------------
# Main
# --------------------------------------------------

if __name__ == "__main__":

    args = parse_args()

    fasta_dir = os.path.join(args.outdir, "fasta")
    protein_dir = os.path.join(fasta_dir, "protein")
    maf_outdir = os.path.join(args.outdir, "maf")

    for d in [fasta_dir, protein_dir, maf_outdir]:
        os.makedirs(d, exist_ok=True)

    species = pd.read_csv(args.species, sep="\t")

    orfs, iorfs, _ = parse_gtf(args.gtf)
    sequences = get_maf(
        iorfs,
        args.maf_dir,
        args.ref_genome,
        fasta_dir,
        maf_outdir
    )

    process_start_codon(sequences, species).to_csv(
        os.path.join(args.outdir, "start_codon_conservation.txt"), sep="\t"
    )
    process_stop_codon(sequences, species).to_csv(
        os.path.join(args.outdir, "stop_codon_conservation.txt"), sep="\t"
    )
    process_sequences_inframe(sequences, species).to_csv(
        os.path.join(args.outdir, "inframe.txt"), sep="\t"
    )
    process_sequences_pid(sequences, species).to_csv(
        os.path.join(args.outdir, "percentage_ids.txt"), sep="\t"
    )

    export_protein_fastas(sequences, protein_dir)
