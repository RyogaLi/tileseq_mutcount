#!/usr/bin/python3.6

import pandas as pd
import numpy as np
import cigar
import os
import sys
import re
import math
import difflib
import itertools
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Seq import MutableSeq


class MutParser(object):

    def __init__(self, row, full_seq, cds_seq, seq_lookup):
        """
        row: input row includes both reads from sam file
        full_seq: full sequence (including cds and padding sequence)
        cds_seq: coding sequence (includes stop codon)
        seq: read from json file - seq, cds_start, cds_end 
        """
        self._seq = Seq(full_seq["seq"].item())
        self._cds = Seq(cds_seq)
        self._start_pos = full_seq.cds_start.item() # 1 posisiton
        self._end_pos = full_seq.cds_end.item()

        # Note all the position in this df are 1-based
        self._seq_lookup = seq_lookup

        # get paired reads and their cigar scores
        self._reads = row

    def _get_seq(self):
        """
        Parse the sequence and CIGAR string
        For R1 and R2, get reference sequence
        """
        # get information for R1
        # where the mapping starts
        self._r1_pos = self._reads.pos_start_r1
        # CIGAR string for r1
        self._r1_cigar = list(cigar.Cigar(self._reads.CIGAR_r1).items())
        # 
        self._r1_readlen = sum([i[0] for i in self._r1_cigar])
        # reference sequence for r1
        self._r1_ref = self._seq[int(self._r1_pos)-1:int(self._r1_pos)+len(self._reads.seq_r1)]
        # quality score string for R1
        self._r1_qual = self._reads.qual_r1
        # R1 sequence
        self._r1_read = self._reads.seq_r1
        # MD:Z string for R1
        self._r1_mdz = self._reads.mdz_r1

        # get the ref sequence for R2
        # where the mapping starts
        self._r2_pos = self._reads.pos_start_r2
        # CIGAR string for R2
        self._r2_cigar = list(cigar.Cigar(self._reads.CIGAR_r2).items())
        #
        self._r2_readlen = sum([i[0] for i in self._r2_cigar])
        # get reference sequence
        self._r2_ref = self._seq[int(self._r2_pos)-1:int(self._r2_pos)+len(self._reads.seq_r2)-1]
        # quality socre string for R2
        self._r2_qual = self._reads.qual_r2
        # R2 sequence
        self._r2_read = self._reads.seq_r2
        # MD:Z string 
        self._r2_mdz = self._reads.mdz_r2

        return self

    def _parse_mut(self, cigar, mdz, ref, read, pos, qual):
        """
        Given a read and information from sam files 
        Return mutations found in this read 
        
        if you are unclear about this function 
        check wiki page: 

        cigar: cigar string
        mdz: mdz string
        ref: Reference sequence for this read
        read: Read sequence
        pos: Starting position of this read
        qual: Quality string

        """
        # join cigar into string
        cigar_joined = [str(i[0])+str(i[1]) for i in cigar]
        pos = int(pos) 
        mut_list = []
        if "I" in "".join(cigar_joined):
            # insertions are only represented in CIGAR string
            # combine mdz and CIGAR to get all the mutations 
            mut_list = self._parse_cigar_ins(cigar, mdz, ref, read, pos, qual)
        else: 
            # when no insertion found, we can process the MD:Z and 
            # find mutations
            mut_list = self._parse_mdz(cigar, mdz, ref, read, pos, qual)
        print("-------------")
        print(mut_list)
        print(cigar)
        print(mdz)
        print(ref)
        print(read)
        print(pos)
        return mut_list

    def _parse_cigar_ins(self, cigar, mdz_raw, ref, read, pos, qual):
        """
        Parse CIGAR string and mdz tags together to identify mutations in a read with insertions
        
        cigar: cigar string
        mdz_raw: mdz string
        ref: Reference sequence for this read
        read: Read sequence
        pos: Starting position of this read
        qual: Quality string
        """
        clip = 0
        # soft clip reads
        # check starting point of the sequence 
        # remove soft clipped bases and adjust for location 
        if cigar[0][1] =="S":
            # adjust read position
            clip += cigar[0][0]

        # convert mdz to list
        mdz = re.findall('.*?[.ATCG]+', mdz_raw)
        
        mut_list = []

        # check insertion poistion in cigar string
        # this postition is corresponded to the reference sequence
        # i.e if pos = 110 means the base 111 is deleted (all 1-based indices)
        ins_pos_ref = 0 # on ref sequence where was the insertion
        total_len = 0 # total lenth indicates in CIGAR string
        ins_pos = [] # store [ins_pos, ins_len]

        mapped = 0
        deleted = 0
        for i in cigar:
            if i[1] != "I": # not insertion
                ins_pos_ref += int(i[0])
                if i[1] == "M": # track number of bases mapped
                    mapped += int(i[0])
                if i[1] == "D": # track number of bases deleted
                    deleted += int(i[0])
            else:
                # based on number of bases mapped, get the inserted base
                # from the read
                ins_base = read[mapped:mapped+int(i[0])]
                # keep the insertion position and inserted lenth in a list
                ins_pos.append([ins_pos_ref, int(i[0])])
                # add the insertion to mut_list
                mut_list.append(str(pos+mapped+deleted)+"|"+ins_base+"|ins")
            total_len += int(i[0])
        
        # parse snp and deletion from mdz string
        # given the inserted position on reference sequence
        r = re.compile("([0-9]+)([a-zA-Z\^]+)")
        
        read_pos = 0 + clip # on the sam read
        deleted_len = 0 # how many bp was deleted
        map_pos = 0 # how many bp was mapped (non insertion track)
        
        iter_ins = iter(ins_pos)
        ins = next(iter_ins, None)
        
        inserted_pos = 0
        for i in mdz:
            # for each item in the MD:Z string 
            # split the item into number and letter
            m = r.match(i)
            match_len = int(m.group(1))
            base = m.group(2)

            map_pos += match_len # update how many bp are mapped
            print(map_pos) 
            
            if ins and map_pos >= ins[0]:
                inserted_pos += ins[1]
                ins = next(iter_ins, None)
            
            if "^" not in base:
                # this means a single nt change
                read_pos += match_len 
                mut_list.append(base+"|"+str(pos+read_pos-clip)+"|"+read[read_pos+inserted_pos-deleted_len])
                map_pos += len(base)
                read_pos += 1
            else: # deletion
                read_pos += match_len
                mut_list.append(str(pos+read_pos-clip)+"|"+base[1:]+"|del")
                deleted_len += len(base[1:])
                read_pos += len(base[1:])
                map_pos -= len(base[1:])
        
        return mut_list

    def _parse_mdz(self, cigar, mdz_raw, ref, read, pos, qual):
        """
        cigar: CIGAR string from SAM file
        mdz: MD:Z:tag from SAM file
        ref: Reference sequence (templete sequence)
        read: Read from SAM file
        read_pos: Mapping position of this read. Where the first base maps to
        qual: Quality score sequence of this read
        
        Get mutation from MD:Z string
        Insertions are not represented by MD:Z tag
        Insertions are processed from CIGAR string
        """
        clip = 0
        # 1. check starting point of the sequence 
        # remove soft clipped bases and adjust for location 
        if cigar[0][1] =="S":
            # adjust read position
            clip += cigar[0][0]
        
        # 2. for the MD:Z string
        # split the string into chunks with mutations
        # i.e ['13C', '14T', '0T', '0G']
        mdz = re.findall('.*?[.ATCG]+', mdz_raw)
        r = re.compile("([0-9]+)([a-zA-Z\^]+)")
        
        read_pos = 0 + clip
        mut_list = []
        deleted_len = 0
        for i in mdz:
            # for each item in the MD:Z string 
            # split the item into number and letter
            m = r.match(i)
            match_len = int(m.group(1))
            base = m.group(2)
            if "^" not in base:
                # this means a single nt change
                # base = reference base
                # pos+read_pos-clip = read starting point + # of bp mapped since that point - clipped base (because they are not on reference sequence)
                # read_pos-deleted_len = # of bases mapped since the beginnig of this read - number of bases deleted
                read_pos += match_len # update read position as we are moving to the right
                # this means a single nt change
                mut_list.append(base+"|"+str(pos+read_pos-clip)+"|"+read[read_pos-deleted_len])
                read_pos += 1 # update read position

            else: # deletion
                read_pos += match_len# update read position as we are moving to the right
                # pos+read_pos-clip = read starting point + # of bp mapped since that point - clipped base (because they are not on reference sequence)
                # base[1:] = bases that were deleted
                mut_list.append(str(pos+read_pos-clip)+"|"+base[1:]+"|del")
                deleted_len += len(base[1:])
                read_pos += len(base[1:])
        
        return mut_list

    def _get_hgvs(self, mut_list):
        """
        mut_list mutation list in the format ['1160|A|1161|del', '1163|T|1164|ins']
        translate a list of mut to aa changes
        return list of aa changes 
        return list of nt changes
        """
        
        # go through mutation list generated by _parse_mdz and _parse_cigar
        # make changes to the DNA sequence 
        # output protein changes for this read
        
        # there are two types of mutations in mut list
        # ins del 
        mut_seq = self._cds.tomutable()
        mutations = ""
        for mut in mut_list:
            if "del" in mut:
                # check reference 
                mut = mut.split("|")
                ref_pos = int(mut[0])
                ref_bp = str(mut[1])
                # get cds position for this mutation 
                cds_pos = self._seq_lookup[self._seq_lookup.temp_pos == ref_pos].cds_pos.item()
                cds_bp = self._seq_lookup[self._seq_lookup.temp_pos == ref_pos].cds.item()
                
                # remove base from mut_seq
                
                
            elif "ins" in mut:
                pass
            else: # snp
                mut_change = mut.split("|")
                hgvs = f"{mut_change[1]}{mut_change[0]}>{mut_change[2]}"
                mutations.append(hgvs)

        return mut_list
    

if __name__ == "__main__":
    ref_dir = "/home/rothlab/rli/02_dev/11_tileSeq/tileseq_py/MTHFR_3del/bowtieIndex/"
    reads = ""
    gene = "MTHFR"

    parser = MutParser(ref_dir, reads, gene)
    parser._parse_mdz([(1,"S"),(140,"M")], "T2A1G", "TGAACG", "CGACCT", 1, "")

    
