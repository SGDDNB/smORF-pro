
smorfs = {} 

#import MySQLdb
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

    with open("smORF_7767.gtf") as gtf:
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

            
def get_maf(iorfs, species):
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
        idx = idxes["chr%s" %(chrom)]
        if iorfs[iorf_id][0][3] == "+": strand = 1
        else: strand = -1
        splices = idx.get_spliced([ int(i[1]) for i in iorfs[iorf_id]], [ int(i[2]) for i in iorfs[iorf_id] ], strand=strand )
        tmp = {}
        for ot in splices: # change to allow dictionary access to a MultipleSeqAlignment
            tmp[ot.id.split(".")[0]] = ot.seq
        all_sequences[iorf_id] = tmp
        
        AlignIO.write(splices, "fasta/%s.fa" %(iorf_id), "fasta")
        AlignIO.write(splices, "maf/%s.maf" %(iorf_id), "maf")
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

def protein_percentageid(sequences, species):
    pid_df = pd.DataFrame(index=list(sequences.keys()), columns=species["GENOME"].tolist())
    for iorf_id in list(sequences.keys()):
        print(iorf_id)
        file1 = open("fasta/protein/%s.txt" %(iorf_id),"a+")
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

import dnds

def dnds(sequences, species):
    dnds_df = pd.DataFrame(index=list(sequences.keys()), columns=species["GENOME"].tolist())
    for iorf_id in list(sequences.keys()):
        print(iorf_id)
        for sp in list(sequences[iorf_id].keys()):
            print(sp)
            seq1 = sequences[iorf_id]["hg38"].replace("-","").translate()
            if len(sequences[iorf_id][sp].replace("-","")) % 3 == 0:
                seq2 = sequences[iorf_id][sp].replace("-","").translate()
            else: continue
            if len(seq1) == len(seq2):
                x = dnds.dnds(str(sequences[iorf_id]["hg38"].replace("-",""))[3:].upper(),
                           str(sequences[iorf_id][sp].replace("-",""))[3:].upper())
                dnds_df.at[iorf_id, sp]= float(x)
    return dnds_df

        
if __name__=="__main__":

    species = pd.read_csv("species_names.txt", sep="\t")
    orfs, iorfs, orf_intervals = parse_gtf("Collapsed_hubner_fixed_V2.gtf")
    sequences = get_maf(iorfs, species)
    
    start_codons = process_start_codon( sequences, species)
    start_codons.to_csv("start_codon_conservation.txt", sep="\t")

    percentage_ids = process_sequences_pid(sequences, species)
    percentage_ids.to_csv("percentage_ids.txt", sep="\t")

    stop_codons = process_stop_codon(sequences, species)
    stop_codons.to_csv("stop_codon_conservation.txt", sep="\t")

    inframe = process_sequences_inframe(sequences, species)
    inframe.to_csv("inframe.txt", sep="\t")

    protein_percentage_ids = protein_percentageid(sequences, species)
    protein_percentage_ids.to_csv("protein_percentage_ids.txt", sep="\t")

    #dnds = bio_parse_orf_gtf.dnds(sequences, species) # need to check if this works ok
    #dnds.to_csv("dnds.txt", sep="\t")
    # amino acid conservation
    
"""
starts = [ int(i[1]) ]
ends = [ int(i[2]) ]
res = idx.search(starts, ends)
fetched = [ multiseq for multiseq in res ]
expected_letters = sum([end - start for start, end in zip(starts, ends)])
all_seqnames = set([sequence.id for multiseq in fetched for sequence in multiseq])
split_by_position = dict([(seq_name, {}) for seq_name in all_seqnames])
positions_test = dict([(seq_name, {}) for seq_name in all_seqnames])
total_rec_length = 0
for multiseq in fetched:
    for seqrec in multiseq:
        if seqrec.id == "hg38.chr8":
            try:
                if ref_first_strand is None:
                    ref_first_strand = seqrec.annotations["strand"]
                    if ref_first_strand not in (1, -1):
                        raise ValueError("Strand must be 1 or -1")
                    elif ref_first_strand != seqrec.annotations["strand"]:
                        raise ValueError("Encountered strand='%s' on target seqname, ""expected '%s'" %(seqrec.annotations["strand"], ref_first_strand))
            except KeyError:
                raise ValueError("No strand information for target seqname (%s)" %self._target_seqname)
            rec_length = len(seqrec)
            rec_start = seqrec.annotations["start"]
            rec_end = seqrec.annotations["start"] + seqrec.annotations["size"]
            total_rec_length += rec_end - rec_start
            for seqrec in multiseq:
                for pos in range(rec_start, rec_end):
                    split_by_position[seqrec.id][pos] = ""
                    positions_test[seqrec.id][pos] = []
            break
        else:
            raise ValueError("Did not find %s in alignment bundle" % (self._target_seqname,))

real_pos = rec_start

for gapped_pos in range(0, rec_length):
    for seqrec in multiseq:
        if seqrec.id == "hg38.chr8":
            print seqrec.seq
            track_val = seqrec.seq[gapped_pos]
            positions_test[seqrec.id][real_pos].append(seqrec.annotations["start"] + (real_pos - seqrec.annotations["start"]))
        print seqrec.id, gapped_pos, track_val, seqrec
        positions_test[seqrec.id][real_pos].append(seqrec.annotations["start"])
        split_by_position[seqrec.id][real_pos] += seqrec.seq[gapped_pos]
    if track_val != "-" and real_pos < rec_end - 1:
        real_pos += 1

        
realpos_to_len = dict([(x, len(y)) for x, y in split_by_position["hg38.chr8"].items() if len(y) > 1])
                                                                                                                                                        
"""
