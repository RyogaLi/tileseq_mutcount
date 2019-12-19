import pandas as pd
import numpy as np
import os
import sys
import glob
import parameters
import pysam 
import parse_seq
import locate_mut
# read from each sam file
# count from pair of sam files

class readSam(object):

    def __init__(self, sam_R1, sam_R2, gene_seq, gene, mut_dir, ref_dir):
        """
        sam_R1: read one of the sample
        sam_R2: read two of the sample
        gene_seq: file contains gene sequence
        gene: gene name 
        mut_dir: directory to save mutation output file
        ref_dir: directory with fasta file
        """
        self._mut_dir = mut_dir
        self._ref_dir = ref_dir
        self._r1 = sam_R1
        self._r2 = sam_R2
        self._sample_id = os.path.basename(self._r1).split("_")[0]
        self._seq_df = gene_seq
        self._gene = gene

    def _trim_sam(self, input_sam):
        """
        trim sam file, convert it to a df with the 
        following colmns:
        read_name, mapped_name, pos_start, mapQ, CIGAR, seq
        """
        trimmed = []
        # read r1 sam file
        with open(input_sam, "r") as sam_file:
            
            for line in sam_file:
                if line.startswith("@"): continue
                line = line.split("\t")
                # create table to save information from sam file
                read_name = line[0]
                mapped_name = line[2]
                pos_start = line[3]
                mapQ = line[4]
                CIGAR = line[5]
                seq = line[9]
                quality = line[10]
                mdz = [i for i in line if "MD:Z:" in i]
                if len(mdz) != 0:
                    mdz = mdz[0].split(":")[-1]
                else:
                    mdz = ""
                trimmed.append([read_name, mapped_name, pos_start,mapQ, CIGAR, mdz, seq, quality])
        
        trimmed = pd.DataFrame(trimmed, columns=["read_name", "mapped_name", "pos_start", "mapQ", "CIGAR", "mdz","seq", "qual"])
        return trimmed

    def _check_qual(self, r1_trimmed, r2_trimmed):
        """
        merged: dataframe that merged R1 and R2 sam files
        """
        # check how many reads mapped in both R1 and R2 (both side have to be high quality of a pair)
        # r1_mapped = r1_trimmed[(r1_trimmed.mapped_name == self._gene) & (r1_trimmed.mapQ.astype(int) >= parameters.mapping_qual)]
        # r2_mapped = r2_trimmed[(r2_trimmed.mapped_name_r2 == self._gene) & (r2_trimmed.mapQ_r2.astype(int) >= parameters.mapping_qual)]
        r1_mapped = r1_trimmed[r1_trimmed.mapped_name == self._gene]
        r2_mapped = r2_trimmed[r2_trimmed.mapped_name_r2 == self._gene]

        # count the number of mapped reads in R1 and R2
        total = r1_trimmed.shape[0]
        r1_n = r1_mapped.shape[0] 
        r2_n = r2_mapped.shape[0]

        # read summary for r1 and r2 are saved to this file
        # sample_idreport.txt
        report_file = self._mut_dir + self._sample_id + "report.txt"
        with open(report_file, "w") as f: 
            f.write("Total number of reads: "+str(total)+"\n")
            f.write("Total number of mapped R1 reads: "+str(r1_n)+"\n")
            f.write("Total number of mapped R2 reads: "+str(r2_n)+"\n")
        
        return r1_mapped, r2_mapped

    def _count(self, mapped_df):
        """
        for each read in mapped_df, count mutations in the reads
        """
        nt_summary = {}
        aa_summary = {}
        
        mut_summary = {"nt_ref":[], "nt_pos":[], "nt_alt":[], "codon_ref":[], "codon_alt":[], "aa_ref":[], "aa_alt":[], "aa_pos":[]}
        mut_df = pd.DataFrame(mut_summary)
        for index, row in mapped_df.iterrows():
            #ps = parse_seq.SeqParser(self._ref_dir, self._gene)
            #ps._map_loc() 
            lm = locate_mut.MutParser(self._ref_dir, row, self._gene)
            read_sum = lm._get_seq()
            
            # compare mutations found in r1 and r2
            # take only the mutations found in both reads
            # alternative implementation: take phred score and 
            # calculate probabilities
            r1_mut_list, r1_mut_seq = lm._parse_mdz(read_sum._r1_cigar, read_sum._r1_mdz, read_sum._r1_ref, read_sum._r1_read, read_sum._r1_pos, read_sum._r1_qual)
            r2_mut_list, r2_mut_seq = lm._parse_mdz(read_sum._r2_cigar, read_sum._r2_mdz, read_sum._r2_ref, read_sum._r2_read, read_sum._r2_pos, read_sum._r2_qual)
            # merge two lists with mutations
            merged_mut = list(set(r1_mut_list)&set(r2_mut_list))
            
            if merged_mut == []:
                continue
            # from merged mutation list, get aa changes
            mutations = lm._translate_mut(merged_mut)
            
            #if mutations == None: 
                # not snp
            #    continue

            mut_df = mut_df.append(mutations)
            #if r1_mut == r2_mut:
            #    print(read_sum._r1_cigar, read_sum._r1_pos)
            #    print(read_sum._r2_cigar, read_sum._r2_pos)
                #print(read_sum._r1_ref)
                #print(read_sum._r2_ref)
            #    mut_list = lm._translate_mut(r1_mut)
            #    print(mut_list)
            #else: # get intersect of two lists
                #mut = tuple(set(r1_mut) & set(r2_mut))
                # print("diff mut")
            #    continue

            #if mut_list:
                #print(mut_list)
            #    for mut in mut_list:
            #        if mut[0] not in nt_summary.keys():
            #            nt_summary[mut[0]]=1
            #        else:
            #            nt_summary[mut[0]]+=1
            #        if mut[1] not in aa_summary.keys():
            #            aa_summary[mut[1]]=1
            #        else:
            #            aa_summary[mut[1]]+=1

        mut_df.to_csv("./mut_df.csv", index=False)
        #break

    def _count_mut(self):
        """
        Starting from sam files (R1 and R2),
        count mutations in each paired reads,
        nt only keep the mutations that are found in both reads
        output files:

        """
        r1_trimmed = self._trim_sam(self._r1)
        r2_trimmed = self._trim_sam(self._r2)
        r2_trimmed.columns = ["read_name", "mapped_name_r2", "pos_start_r2", "mapQ_r2", "CIGAR_r2", "mdz_r2", "seq_r2", "qual_r2"]
        
        # check quality of the reads
        # mapping criteria: the mapped_name must equal to gene name
        # mapQ greater than quality threshold (see parameters.py)
        r1_mapped, r2_mapped = self._check_qual(r1_trimmed, r2_trimmed)
        
        # merge R1 and R2 reads into one df
        # for each pair of R1 and R2 reads, we need to see the mutation
        # on both sequences
        merged = pd.merge(r1_mapped, r2_mapped, on="read_name")
        
        # from merged, filtered reads, count mutations 
        self._count(merged)

