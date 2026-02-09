#!/usr/bin/env python3
"""
Conservation analysis script for smORF-pro
Modified to accept command-line arguments
"""

import argparse
import os
import sys

smorfs = {} 

# import MySQLdb
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq

from bx.align import maf
import pandas as pd
import subprocess
import re
import functools


def parse_gtf(fn):
    orfs = {}
    iorfs = {}
    orf_intervals = {}

    with open(fn) as gtf:
        for l in gtf:
            fields = l.split("\t")
            interval = ( fields[0], str(int(fields[3])-1), fields[4], fields[6] )
            attributes = fields[8]
            attributes_split = dict([ i.split(" ") for i in attributes.split("; ") ])
            if fields[2] == "iORF":
                orfs[ attributes_split["ORF_id"] ] = set()
                orf_intervals[ attributes_split["ORF_id"] ] = interval 
            elif fields[2] == "orfCDS":
                orfs[ attributes_split["ORF_id"] ].add(  attributes_split["iORF_id" ] )
                if attributes_split["iORF_id" ] not in list(iorfs.keys()):
                    iorfs[ attributes_split["iORF_id" ] ] = []
                iorfs[ attributes_split["iORF_id" ] ].append( interval )
    print(len(list(orfs.keys())))
    print(len(list(iorfs.keys())))
    for iorf_id in list(iorfs.keys()):
        iorfs[ iorf_id ].sort(key=functools.cmp_to_key(sort_intervals))
    return orfs, iorfs, orf_intervals
    
def sort_intervals(x, y):
    if x[3] == "+":
        return 1 if x[1] > y[1] else 0 if x[1] == y[1] else -1
    else:
        return 1 if x[2] > y[2] else 0 if x[2] == y[2] else -1


def get_locations(iorfs, species):
    """ this was the start of attempting to extract the start and end locations of aligned smorfs"""
    all_sequences = {}
    idxes = {}
    for chr in [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                 "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                 "chrX", "chrY", "chrM" ]:
        idxes[chr] = AlignIO.MafIO.MafIndex("%s.mafindex" %(chr), "/home/baptiste/Sonia_computer/%s.maf" %(chr), "hg38.%s" %(chr))
        for iorf_id in list(iorfs.keys()):
            iorf = iorfs[iorf_id]
            chrom = iorf[0][0]
            region_name = "%s.chr%s" % ("hg38", chrom)
        else: strand = -1
        splices = idx.get_spliced([ int(i[1]) for i in iorfs[iorf_id]], [ int(i[2]) for i in iorfs[iorf_id] ], strand=strand )
        tmp = {}
        for ot in splices: # change to allow dictionary access to a MultipleSeqAlignment
            tmp[ot.id.split(".")[0]] = ot.seq
            all_sequences[iorf_id] = tmp
            
        """AlignIO.write(splices, "fasta/%s.fa" %(iorf_id), "fasta")"""
    return all_sequences

            
def get_maf(iorfs, species, maf_dir, fasta_dir, protein_dir):
    all_sequences = {} 
    idxes = {}
    
    # Get list of available chromosomes from MAF directory
    available_chr = []
    for chr in [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
               "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
               "chrX", "chrY", "chrM" ]:
        maf_file = os.path.join(maf_dir, "%s.maf" % chr)
        if os.path.exists(maf_file):
            available_chr.append(chr)
            idx_file = os.path.join(maf_dir, "%s.mafindex" % chr)
            idxes[chr] = AlignIO.MafIO.MafIndex(idx_file, maf_file, "hg38.%s" % chr)
    
    print("Available chromosomes:", available_chr)
    
    for iorf_id in list(iorfs.keys()):
        iorf = iorfs[iorf_id]
        chrom = iorf[0][0]
        # Normalize chromosome name (add 'chr' prefix if missing)
        chr_key = chrom if chrom.startswith("chr") else "chr%s" % chrom
        
        if chr_key not in idxes:
            print("Warning: MAF file not found for chromosome %s (from %s), skipping iORF %s" % (chr_key, chrom, iorf_id))
            continue
            
        region_name = "%s.%s" % ("hg38", chr_key)
        idx = idxes[chr_key]
        if iorfs[iorf_id][0][3] == "+": strand = 1
        else: strand = -1
        
        try:
            splices = idx.get_spliced([ int(i[1]) for i in iorfs[iorf_id]], [ int(i[2]) for i in iorfs[iorf_id] ], strand=strand )
            tmp = {}
            for ot in splices: # change to allow dictionary access to a MultipleSeqAlignment
                tmp[ot.id.split(".")[0]] = ot.seq
            all_sequences[iorf_id] = tmp
            
            AlignIO.write(splices, os.path.join(fasta_dir, "%s.fa" % iorf_id), "fasta")
        except Exception as e:
            print("Error processing iORF %s: %s" % (iorf_id, str(e)))
            continue
    
    return all_sequences

# sequences[iorf_id]["hg38"].ungap("-").translate()

def process_start_codon(sequences, species):
    sc_df = pd.DataFrame(index=list(sequences.keys()), columns=species["GENOME"].tolist())
    for iorf_id in list(sequences.keys()):
        seqs = [ x.start() for x in re.finditer("[atgcn]", str(sequences[iorf_id]["hg38"]).lower() ) ]
        first_codon_indexes = seqs[0:3]
        seq = str(sequences[iorf_id]["hg38"])[:(first_codon_indexes[2]+1)]
        gaps = [ x.start() for x in re.finditer("-", seq) ]
        human_start_codon = "".join([str(sequences[iorf_id]["hg38"]).lower()[index] for index in first_codon_indexes ])
        for sp in list(sequences[iorf_id].keys()):
            tmp_seq = str(sequences[iorf_id][sp])[:(first_codon_indexes[2]+1)]
            tmp_gaps = [ x.start() for x in re.finditer("-", tmp_seq) ]
            if gaps == tmp_gaps:
                start_codon = "".join([str(sequences[iorf_id][sp]).lower()[index] for index in first_codon_indexes ])
            else:
                print("HERE", gaps, tmp_gaps, sp, iorf_id)
                start_codon = str(sequences[iorf_id][sp]).lower()[:3]
            #sc_df.set_value(iorf_id, sp, start_codon == "atg")
        sc_df.at[iorf_id,sp]=start_codon
    return sc_df

def process_stop_codon(sequences, species):
    sc_df = pd.DataFrame(index=list(sequences.keys()), columns=species["GENOME"].tolist())
    for iorf_id in list(sequences.keys()):
        print(iorf_id)
        seqs = [ x.start() for x in re.finditer("[atgcn]", str(sequences[iorf_id]["hg38"]).lower() ) ]
        last_codon_indexes = seqs[(len(seqs)-3):len(seqs)]
        human_stop_codon = "".join([str(sequences[iorf_id]["hg38"]).lower()[index] for index in last_codon_indexes ])
        seq = str(sequences[iorf_id]["hg38"])[last_codon_indexes[0]:]
        gaps = [ x.start() for x in re.finditer("-", seq) ]
        for sp in list(sequences[iorf_id].keys()):
            tmp_seq = str(sequences[iorf_id][sp])[last_codon_indexes[0]:]
            tmp_gaps = [ x.start() for x in re.finditer("-", tmp_seq) ]
            if gaps == tmp_gaps:
                 stop_codon = "".join([str(sequences[iorf_id][sp]).lower()[index] for index in last_codon_indexes ])
            else:
                stop_codon = str(sequences[iorf_id][sp]).lower()[len(str(sequences[iorf_id][sp]))-3:]
            sc_df.at[iorf_id, sp]=stop_codon
    return sc_df
        
def process_sequences_inframe(sequences, species):
    if_df = pd.DataFrame(index=list(sequences.keys()), columns=species["GENOME"].tolist())
    for iorf_id in list(sequences.keys()):
        for sp in list(sequences[iorf_id].keys()):
            if len(sequences[iorf_id][sp].replace("-","")) % 3 == 0:   
                if_df.at[iorf_id, sp]=True
            else:
                if_df.at[iorf_id, sp]=False
    return if_df

def process_sequences_pid(sequences, species, target="hg38"):
    pid_df = pd.DataFrame(index=list(sequences.keys()), columns=species["GENOME"].tolist())
    for iorf_id in list(sequences.keys()):
        print(iorf_id)
        s = [ s for s in species["GENOME"].tolist() if s in list(sequences[iorf_id].keys()) ]
        print(s)
        for i in s:
            print(i)
            if sequences[iorf_id][i].count("-") == len(sequences[iorf_id][target]):
                pid_df.at[iorf_id, i]=None
            # TODO the following condition is not required anymore!
            if len(sequences[iorf_id][target]) == len(sequences[iorf_id][i]) and sequences[iorf_id][i].count("-") != len(sequences[iorf_id][i]):
                match = 0
                total = 0
                for ic in range( 0, len(sequences[iorf_id][target]) ):
                    a = sequences[iorf_id][target][ic].lower()
                    b = sequences[iorf_id][i][ic].lower()
                    if a == '-' or b == '-': # gaps are not a match
                        continue
                    elif a == b:
                        match += 1
                    total += 1
                if total == 0:
                    pass
                else: pid_df.at[iorf_id, i]=100 * ( match / float(total))
    return pid_df

def protein_percentageid(sequences, species, protein_dir):
    pid_df = pd.DataFrame(index=list(sequences.keys()), columns=species["GENOME"].tolist())
    for iorf_id in list(sequences.keys()):
        print(iorf_id)
        file1 = open(os.path.join(protein_dir, "%s.txt" % iorf_id), "w")
        #file1.write("\n>hg38\n")
        #for listitem in seq1:
        #        file1.write('%s' % listitem)
        for sp in list(sequences[iorf_id].keys()):
            print(sp)
            seq1 = sequences[iorf_id]["hg38"].replace("-","").translate()
            if len(sequences[iorf_id][sp].replace("-","")) % 3 == 0:
                seq2 = sequences[iorf_id][sp].replace("-","").translate()
            else: continue

            match = 0
            total = 0
            #for ic in range( 0, len(seq1) ):
            #    a = seq1[ic].lower()
            #    b = seq2[ic].lower()
            #    if a == '-' or b == '-': # gaps are not a match
            #        continue
            #    elif a == b:
            #        match += 1
            #    total += 1
            #if total == 0:
            #    pass
            #else: 
            #    pid_df.at[iorf_id, sp]=100 * ( match / float(total))
            file1.write('\n>%s\n' % sp)
            for listitem in seq2:
                file1.write('%s' % listitem)
        file1.close()
    return pid_df

# Note: dnds function commented out as it requires additional package
# import dnds
# def dnds(sequences, species):
#     ...

        
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Conservation analysis for smORF-pro')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--maf-dir', required=True, help='Directory containing MAF files')
    parser.add_argument('--output-dir', required=True, help='Output directory (will create fasta and fasta/protein subdirectories)')
    parser.add_argument('--species', required=True, help='Path to species_names.txt file')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.gtf):
        print("Error: GTF file not found: %s" % args.gtf)
        sys.exit(1)
    
    if not os.path.isdir(args.maf_dir):
        print("Error: MAF directory not found: %s" % args.maf_dir)
        sys.exit(1)
    
    if not os.path.exists(args.species):
        print("Error: Species names file not found: %s" % args.species)
        sys.exit(1)
    
    # Create output directories
    fasta_dir = os.path.join(args.output_dir, "fasta")
    protein_dir = os.path.join(args.output_dir, "fasta", "protein")
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(protein_dir, exist_ok=True)
    
    # Change to output directory for file writing
    original_cwd = os.getcwd()
    os.chdir(args.output_dir)
    
    try:
        # Load species names
        species = pd.read_csv(args.species, sep="\t")
        
        # Parse GTF
        print("Parsing GTF file: %s" % args.gtf)
        orfs, iorfs, orf_intervals = parse_gtf(args.gtf)
        
        # Get MAF alignments
        print("Processing MAF files from: %s" % args.maf_dir)
        sequences = get_maf(iorfs, species, args.maf_dir, fasta_dir, protein_dir)
        
        # Process results
        print("Processing start codons...")
        start_codons = process_start_codon(sequences, species)
        start_codons.to_csv("start_codon_conservation.txt", sep="\t")
        
        print("Processing percentage identity...")
        percentage_ids = process_sequences_pid(sequences, species)
        percentage_ids.to_csv("percentage_ids.txt", sep="\t")
        
        print("Processing stop codons...")
        stop_codons = process_stop_codon(sequences, species)
        stop_codons.to_csv("stop_codon_conservation.txt", sep="\t")
        
        print("Processing inframe sequences...")
        inframe = process_sequences_inframe(sequences, species)
        inframe.to_csv("inframe.txt", sep="\t")
        
        print("Processing protein percentage identity...")
        protein_percentage_ids = protein_percentageid(sequences, species, protein_dir)
        protein_percentage_ids.to_csv("protein_percentage_ids.txt", sep="\t")
        
        print("Conservation analysis complete!")
        
    finally:
        # Restore original working directory
        os.chdir(original_cwd)
