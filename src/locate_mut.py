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
        self._r1_ref = self._seq[int(self._r1_pos)-1:int(self._r1_pos)+len(self._reads.seq_r1)-1]
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
        
        mut_list = []
        if "I" in "".join(cigar_joined):
            # insertions are only represented in CIGAR string
            # combine mdz and CIGAR to get all the mutations 
            print(cigar)
            print(mdz)
            exit()
        else: 
            # when no insertion found, we can process the MD:Z and 
            # find mutations
            mut_list = self._parse_mdz(cigar, mdz, ref, read, pos, qual)
        
        return mut_list

        

    def _parse_mdz(self, cigar, mdz_raw, ref, read, pos, qual):
        """
        cigar: CIGAR string from SAM file
        mdz: MD:Z:tag from SAM file
        ref: Reference sequence (templete sequence)
        read: Read from SAM file
        read_pos: Mapping position of this read
        qual: Quality score sequence of this read
        
        Get mutation from MD:Z string
        Insertions are not represented by MD:Z tag
        Insertions are processed from CIGAR string
        """
        pos = int(pos)

        # 1. check starting point of the sequence 
        # remove soft clipped bases and adjust for location 
        if cigar[0][1] =="S":
            # adjust read position
            pos -= cigar[0][0]
        
        # 2. for the MD:Z string
        # split the string into chunks with mutations
        # i.e ['13C', '14T', '0T', '0G']
        mdz = re.findall('.*?[.ATCG]+', mdz_raw)
        r = re.compile("([0-9]+)([a-zA-Z\^]+)")
        
        read_pos = 0
        mut_list = []
        for i in mdz:
            # for each item in the MD:Z string 
            # split the item into number and letter
            m = r.match(i)
            match_len = int(m.group(1))
            base = m.group(2)
            if "^" not in base:
                # this mean a single nt change
                mut_list.append(base+"|"+str(pos+read_pos+match_len)+"|"+read[read_pos+match_len])
                read_pos += match_len+1 
            else: # deletion
                mut_list.append(str(pos+match_len+1)+"|"+base[1:]+"|del")
                read_pos += match_len
        
        return mut_list

    def _parse_cigar(self, cigar, ref, read, read_pos, qual):
        """
        given a list of cigar string, read and reference read
        find the mutations based on CIGAR string
        cigar: list of cigar strings 
        ref: reference sequence
        read: read in sam file
        read_pos: start position of the read
        qual: quality of the read
        """
        # pos = start position of the read on ref
        start_pos = int(read_pos)-1

        read_inx = 0
        # if the read is soft clipped 
        if cigar[0][1] =="S":
            # adjust read position
            start_pos -= cigar[0][0]
        #if len(cigar) == 1 and cigar[0][1] == "M":
            # the read mapped, no mutation found

        ref_pos = int(start_pos)
        mut_seq = self._seq[:ref_pos]
        mut = []
        n_del = 0
        for i in cigar:
            
            if i[1] == "H":
                # need an example of this read
                ref_pos+=i[0]
                
            if i[1] == "S": # position clipped
                # skipp this region  
                read_inx += i[0]
                ref_pos+=i[0]
            
            if i[1] == "M": # mapped region
                
                read_inx += i[0]
                ref_pos += i[0]

            if i[1] == "D": # deletion 
                # ref base 
                # get the deleted base from reference 
                # store mut in the format: start pos | ref bases | end_pos | type
                ref_bp=self._seq[ref_pos:ref_pos+i[0]]
                nt_mut = str(ref_pos)+"|"+str(ref_bp)+"|"+str(ref_pos+int(i[0]))+"|del"
                ref_pos += i[0]
                #read_inx += i[0]
                mut.append(nt_mut)
                n_del +=1

            if i[1] == "I": # insertion
                
                # from the read, get the bases that are inserted 
                # store the mut in the format: start pos | inserted bases | end pos | type
                inserted_bp = read[read_inx:read_inx+i[0]]
                #print("merged insert")
                mut_seq += inserted_bp
                #print(mut_seq)
                #print(inserted_bp)
                nt_mut = str(ref_pos)+"|"+str(inserted_bp)+"|"+str(ref_pos+i[0])+"|ins"
                ref_pos += len(inserted_bp)-1
                #print(ref_pos)
                mut.append(nt_mut)
                n_del -=1

        # add the rest seq to mut_seq
        # mut_seq += self._seq[start_pos+len(read)+n_del:]
        return mut

    def _convert_qual(self):
        """
        get r1 qual and r2 qual
        """
        pass
    
    def _get_aa(self, pos, ref, alt, t=None):
        """
        pos: nt position
        ref: reference bases 
        alt: alt bases
        """
        if t == "SNP":
            # convert pos to protein loc
            pro_loc = math.ceil(pos / 3)
            mut = self._ref_pro[pro_loc-1]+"|"+str(pro_loc)+"|"+self._mut_pro[pro_loc-1]
        elif t == "del":
            # in case of deletion
            # pos is a tuple contains the nt range
            # alt is NONE

            pro_start = math.ceil(pos[0] / 3)
            pro_end = math.ceil(pos[1] / 3)
            deleted_pro = Seq(ref).translate()

        elif t == "ins":
            # incase of insertion
            # pos is a tuple contains the nt range for the insertion
            # ref is NONE
            pass
        return mut

    def _translate_mut(self, mut_list):
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
            else:
                pass

        return mut_list
    
    def _indel_change(self, ref, mut):
        """
        small insertion and deletions
        """
        pass

    def _mut_change(self, ref, mut):
        """
        len ref and len mut are the same
        given reference seq and mut seq, find out nt changes 
        return nt changes in the format of ref, pos, alt, codon_ref, codon_alt
        """
        
        mut_summary = {"nt_ref":[], "nt_pos":[], "nt_alt":[], "codon_ref":[], "codon_alt":[], "aa_ref":[], "aa_alt":[], "aa_pos":[]}
        nt = 0
        while nt < len(ref):
            if ref[nt] == mut[nt]:
                nt+=1
                continue 
            else:
                # get 3 bp codon of this nt
                if (nt+1) %3 == 0: # last bp of a codon
                    c1_ref, c1_alt = ref[nt-2], mut[nt-2]
                    c2_ref, c2_alt = ref[nt-1], mut[nt-1]
                    c3_ref, c3_alt = ref[nt], mut[nt]
                elif (nt+1)%3 ==2: # second bp of a codon
                    c1_ref, c1_alt = ref[nt-1], mut[nt-1]
                    c2_ref, c2_alt = ref[nt], mut[nt]
                    c3_ref, c3_alt =ref[nt+1], mut[nt+1]
                elif (nt+1)%3 ==1: # first bp of a codon
                    c1_ref, c1_alt = ref[nt], mut[nt]
                    c2_ref, c2_alt = ref[nt+1], mut[nt+1]
                    c3_ref, c3_alt = ref[nt+2], mut[nt+2]
                codon_ref = c1_ref+c2_ref+c3_ref
                codon_alt = c1_alt+c2_alt+c3_alt
                aa_pos = math.ceil((nt+1)/3)
                mut_summary["codon_ref"].append(codon_ref)
                mut_summary["codon_alt"].append(codon_alt)
                mut_summary["aa_ref"].append(str(Seq(codon_ref).translate()))
                mut_summary["aa_alt"].append(str(Seq(codon_alt).translate()))
                mut_summary["aa_pos"].append(aa_pos)
                mut_summary["nt_ref"].append(ref[nt])
                mut_summary["nt_alt"].append(mut[nt])
                mut_summary["nt_pos"].append(nt+1) # 1 based position
            nt+=1 
        return pd.DataFrame(mut_summary)
    
    def _mut_change_mod(self, ref, mut):
        """
        get nt changes
        """
        mut_summary = {"nt_ref":[], "nt_pos":[], "nt_alt":[]}

        nt = 0
        while nt < len(ref):
            if ref[nt] == mut[nt]:
                nt+=1
                continue
            else:
                if (nt+1) %3 == 0: # last bp of a codon
                    codon_ref = ref[nt-2]+ref[nt-1]+ref[nt]
                elif (nt+1)%3 ==2: # second bp of a codon
                    codon_ref = ref[nt-1]+ref[nt]+ref[nt+1]
                elif (nt+1)%3 ==1: # first bp of a codon
                    codon_ref = ref[nt]+ref[nt+1]+ref[nt+2]
                            
            mut_summary["nt_ref"].append(ref[nt])
            mut_summary["nt_alt"].append(mut[nt])
            mut_summary["nt_pos"].append(nt+1)
            mut_summary["aa_ref"].append(str(Seq(codon_ref).translate()))
            

if __name__ == "__main__":
    ref_dir = "/home/rothlab/rli/02_dev/11_tileSeq/tileseq_py/MTHFR_3del/bowtieIndex/"
    reads = ""
    gene = "MTHFR"

    parser = MutParser(ref_dir, reads, gene)
    parser._parse_mdz([(1,"S"),(140,"M")], "T2A1G", "TGAACG", "CGACCT", 1, "")

    
