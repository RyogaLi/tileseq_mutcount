#!/usr/bin/env python3.6

#  What does this script do?
# This script takes a pair of read (R1 and R2) in the format of dataframe row
# and call mutaitons based on the reads


# TODO

import pandas as pd
import cigar
import re
from Bio.Seq import Seq

import posterior

class MutParser(object):

    def __init__(self, row, full_seq, cds_seq, seq_lookup, tile_s, tile_e, post_prob_cutoff, logging, mut_rate):
        """
        row: input row includes both reads from sam file
        full_seq: full sequence (including cds and padding sequence)
        cds_seq: coding sequence (includes stop codon)
        seq: read from json file - seq, cds_start, cds_end
        """
        self._seq = Seq(full_seq["seq"].values.item())
        self._cds = Seq(cds_seq)
        self._start_pos = full_seq.cds_start.values.item() # 1 posisiton
        self._end_pos = full_seq.cds_end.values.item()

        self._tile_begins = tile_s
        self._tile_ends = tile_e
        # Note all the position in this df are 1-based
        self._seq_lookup = seq_lookup

        # get paired reads and their cigar scores
        self._reads = row

        # user defined post_prob_cutoff
        self._cutoff = post_prob_cutoff
        self._mutrate = mut_rate

        # logger
        self._logging = logging



    def _get_seq(self):
        """
        Parse the sequence and CIGAR string
        For R1 and R2, get reference sequence
        """
        # get information for R1
        # where the mapping starts
        self._r1_pos = int(self._reads["pos_start_r1"])
        self._r1_cigar = list(cigar.Cigar(self._reads["cigar_r1"]).items())
        self._r1_readlen = sum([i[0] for i in self._r1_cigar])
        self._r1_ref = self._seq[int(self._r1_pos)-1:int(self._r1_pos)+len(self._reads["seq_r1"])]
        self._r1_qual = self._reads["qual_r1"]
        self._r1_read = self._reads["seq_r1"]
        self._r1_mdz = self._reads["mdz_r1"]

        # get the ref sequence for R2
        self._r2_pos = int(self._reads["pos_start_r2"])
        self._r2_cigar = list(cigar.Cigar(self._reads["cigar_r2"]).items())
        self._r2_readlen = sum([i[0] for i in self._r2_cigar])
        self._r2_ref = self._seq[int(self._r2_pos)-1:int(self._r2_pos)+len(self._reads["seq_r2"])]
        self._r2_qual = self._reads["qual_r2"]
        self._r2_read = self._reads["seq_r2"]
        self._r2_mdz = self._reads["mdz_r2"]

        return self

    def _main(self):
        """
        return a list of mutations from paired reads (R1 and R2)
        """
        # assign names to items in the dictionary
        self._get_seq()

        # parse mutations in R1
        r1_mut = self._parse_cigar_mdz(self._r1_cigar, self._r1_mdz, self._r1_ref, self._r1_read, self._r1_pos, self._r1_qual)
        # parse mutations in R2
        r2_mut = self._parse_cigar_mdz(self._r2_cigar, self._r2_mdz, self._r2_ref, self._r2_read, self._r2_pos, self._r2_qual)

        # for deletion or insertion, we count those that appeared on both reads
        # del_ins = list(set(r1_delins) & set(r2_delins))
        if r1_mut != [] or r2_mut != []:
            final_mut = []
        #if len(r1_mut)>5:
            r1_m = pd.DataFrame({"m_r1": r1_mut})
            r1_m["read"] = "r1"
            r2_m = pd.DataFrame({"m_r2": r2_mut})
            r2_m["read"] = "r2"
            if not r1_mut == []: 
                r1_m = pd.DataFrame({"m_r1": r1_mut})
                r1_m["read"] = "r1"
                r2_m = pd.DataFrame({"m_r2": r2_mut})
                r2_m["read"] = "r2"
                r1_m[["pos", "ref_r1", "alt_r1", "qual_r1"]] = r1_m["m_r1"].str.split("|", expand=True)
            else:
                r1_m = pd.DataFrame([], columns = ["m_r1", "pos", "ref_r1", "alt_r1", "qual_r1", "read"])

            if not r2_mut == []:
                r2_m = pd.DataFrame({"m_r2": r2_mut})
                r2_m["read"] = "r2"
                r2_m[["pos", "ref_r2", "alt_r2", "qual_r2"]] = r2_m["m_r2"].str.split("|", expand=True)
            else: 
                r2_m = pd.DataFrame([], columns = ["m_r2", "pos", "ref_r2", "alt_r2", "qual_r2", "read"])

            snp_df = pd.merge(r1_m, r2_m, on=["pos"], how="outer")

            snp_df["pos"] = pd.to_numeric(snp_df["pos"])
            # group mutations based on positions
            n = 3 #tmp
            d = dict(tuple(snp_df.groupby(snp_df['pos'].diff().gt(n).cumsum())))
            
            pos_df = posterior.cluster(d, self._mutrate, self._cutoff) # analyze the dictionary of clusters and get posterior
            # two cases of snp
            # mutations are on both reads
            # mutaiton only shows on one read

            # d = dict(tuple(df.groupby(df['x'].diff().gt(100).cumsum())))
            # 1. group mutations by position
            # split snp column into 3
            #print(snp_df)
                #pos_prob = posterior.bayesian_variant_call([r1_basecall, r2_basecall], [r1_qual, r2_qual], wt, self._mutrate)
                #print(pos_prob)
                # if two mut are different, it will return a dictionary with two keys and their posterior probabilities
                # pick the one with higher probability
                # if it is the wt then do nothing
                # if not, add this mutation to the final list
                #if pos_prob[r1_basecall] > pos_prob[r2_basecall]:
                #    if r1_basecall == wt: continue
                #    if pos_prob[r1_basecall] > self._cutoff:
                #        final_mut.append(row.snp)
                #elif pos_prob[r1_basecall] < pos_prob[r2_basecall]:
                #    if r2_basecall == wt: continue
                #    if pos_prob[r2_basecall] > self._cutoff:
                #        final_mut.append(row.snp)
                #elif pos_prob[r1_basecall] == pos_prob[r2_basecall]:
                #    if r2_basecall == wt or r1_basecall == wt:
                #        print(r1_basecall, r2_basecall, wt)
                #    if pos_prob[r2_basecall] > self._cutoff:
                #        final_mut.append(row.snp)
                
            final_mut = list(set(pos_df.m.tolist()))
            final_mut.sort()
            if final_mut != []:
                hgvs, outside_mut = self._get_hgvs(final_mut)
            else:
                hgvs, outside_mut, pos_df = [],[],[]
        return hgvs, outside_mut, pos_df

    def _parse_cigar_mdz(self, cigar, mdz_raw, ref, read, pos, qual):
        """
        use CIGAR string and mdz tag to call mutations from a read
        """
        snp_list = [] # output list with all the mutations
        delins_list = []
        # soft clip
        clip = 0
        # check starting point of the sequence
        # remove soft clipped bases and adjust for position
        if cigar[0][1] == "S": # soft clip occurs in the beginning of a read
            clip += cigar[0][0]

        # convert mdz string into a list
        mdz = re.findall('.*?[.ATCG]+', mdz_raw)
        # convert cigar list to a string
        cigar_joined = [str(i[0])+str(i[1]) for i in cigar]

        ins_pos_ref = 0 # on ref sequence where was the insertion
        total_len = 0 # total lenth indicates in CIGAR string
        ins_pos = [] # store [ins_pos, ins_len]

        mapped = 0
        deleted = 0

        # create a dictionary to store the mapped positions (ref_pos:read_pos)
        # this is later used to get the quality score for a base
        map_pos = {}
        ref_positions = []
        read_positions = []
        read_start = 0+clip
        ref_start = pos
        #if "I" in "".join(cigar_joined): # if insertion is in the read
            # get insertion position from cigar string
        for i in cigar: # go through each cigar string
            if i[1] != "I": # not insertion
                ins_pos_ref += int(i[0])
                if i[1] == "M": # track number of bases mapped
                    ref_positions += list(range(ref_start, ref_start+int(i[0])))
                    read_positions += list(range(read_start, read_start+int(i[0])))
                    ref_start += int(i[0])
                    read_start += int(i[0])
                    mapped += int(i[0])

                if i[1] == "D": # track number of bases deleted
                    deleted += int(i[0])
                    #read_start += int(i[0])
                    # skip deleted bases on ref read
                    ref_start += int(i[0])

            else:
                # based on number of bases mapped, get the inserted base from read
                ins_base = read[read_start:read_start+int(i[0])]
                qual_base = qual[read_start:read_start+int(i[0])]
                qual_base = list(qual_base)
                qual_base = ",".join(qual_base)
                # calculate inserted base posterior by using the quality string
                # qual_ins = qual[mapped:mapped+int(i[0])]
                # keep the insertion position and inserted lenth in a list
                ins_pos.append([ins_pos_ref, int(i[0])])
                # add the insertion to mut_list
                #mut_list.append([str(pos+mapped+deleted),ins_base,"ins"])
                delins_list.append(f"{pos+mapped+deleted}|{ins_base}|ins|{qual_base}")
                # skip inserted bases on read
                read_start += int(i[0])
        # zip two lists into a dictionary
        pos_map = dict(zip(ref_positions, read_positions))
        # parse snp and deletion from mdz string
        # given the inserted position on reference sequence
        r = re.compile("([0-9]+)([a-zA-Z\^]+)")
        # create iterator to go through all the insertions
        iter_ins = iter(ins_pos)
        ins = next(iter_ins, None)

        read_pos = 0
        inserted_pos = 0
        deleted_len = 0
        map_pos = 0 + clip
        ##read_ref_map = [] # map [ref_pos, ref_base, read_pos, read_base, read_quality_score]
        for i in mdz:
            # for each item in the MD:Z string
            # split the item into number and letter
            m = r.match(i)
            match_len = int(m.group(1))
            base = m.group(2)

            map_pos += match_len# update how many bp are mapped
            read_pos += match_len

            while ins and map_pos >= ins[0]:
                map_pos += ins[1]
                inserted_pos += ins[1]
                ins = next(iter_ins, None)

            if "^" not in base:
                # this means a single nt change
                #mut_list.append([str(pos+read_pos-clip),base,read[read_pos+inserted_pos-deleted_len],qual[read_pos+inserted_pos-deleted_len]])
                snp_list.append(f"{pos+read_pos}|{base}|{read[read_pos+inserted_pos-deleted_len+clip]}|{qual[read_pos+inserted_pos-deleted_len+clip]}")

                map_pos += len(base)
                read_pos += len(base) # adjust read pos with 1bp change (move to the right for 1 pos)
            else: # deletion
                #mut_list.append([str(pos+read_pos-clip),base[1:],"del"])
                delins_list.append(f"{pos+read_pos}|{base[1:]}|del|{qual[read_pos+inserted_pos-deleted_len+clip-1]},{qual[read_pos+inserted_pos-deleted_len+clip]}")

                deleted_len += len(base[1:])
                read_pos += len(base[1:])
                map_pos += len(base[1:])
        
        if len(delins_list) >= 3:
            print(snp_list, delins_list, pos_map)
            print(mdz)
            print(cigar)
            print(ref)
            print(read)
            print(qual)
            print("---")
            #exit()
        mut_list = snp_list+delins_list
        return mut_list

    def _get_hgvs(self, mut_list):
        """
        mut_list mutation list in the format
        translate a list of mut to aa changes
        return list of nt changes represent by hgvs strings
        """
        # there are two types of mutations in mut list
        # ins del

        # how to track consecutive changes?
        concecutive_snp = [] # if the snp changes are concecutive, it will be represent as delins
        combined_snp = ""

        # track mutation positions that was not within the tiled region
        outside_mut = []

        # use to track delins
        # if del and ins happened within 2bp from each other then we consider it a delins
        # also record any snp that happened in the range
        delins = []
        mutations = []
        for mut in mut_list:
            mut_change = mut.split("|")
            tmp_pos = int(mut_change[0]) # template position

            # get cds position from lookup table
            try:
                cds_pos = self._seq_lookup[self._seq_lookup.temp_pos == tmp_pos].cds_pos.values.item()
            except:
                continue

            if cds_pos < self._tile_begins or cds_pos > self._tile_ends:
                #self._logging.warning(f"mutation at pos {cds_pos} which is not within the tile")
                outside_mut.append(cds_pos)
                continue

            if "N" in mut: # do not consider base N
                continue

            if ("del" in mut) or ("ins" in mut): # deletion or insertion
                mut_change[0] = cds_pos
                if delins == []: # nothing is on track
                    delins.append(mut_change)
                else: # compare bases
                    # delins = [[19,T,del]]
                    prev_pos = delins[-1][0] # previous deletion or insertion position
                    prev_bases = delins[-1][1]
                    prev_change = delins[-1][2]
                    if "del" in mut:
                        record_pos = mut_change[0]
                    elif "ins" in mut:
                        record_pos = mut_change[0] - 1

                    if (prev_pos-2 <= record_pos <= prev_pos+2) and (mut_change[2] != prev_change): # within 2bp of previous change
                        # not concecutive del or ins
                        # add this to delins
                        delins.append(mut_change)
                    else:
                        # output hgvs of mutations in delins
                        hgvs = delins_to_hgvs(self._cds, delins)
                        # update delins with current mutation
                        delins = [mut_change]
                        mutations.append(hgvs)

            else: # snp
                #mut_change = mut.split("|")
                #tmp_pos = int(mut_change[0])
                # get cds position from look up table
                #cds_pos = self._seq_lookup[self._seq_lookup.temp_pos == tmp_pos].cds_pos.item()
                cds_ref = self._cds[cds_pos-1]

                ## validate: if the reference base matches the ref in seq_lookup
                if cds_ref != mut_change[1]:
                    raise ValueError(f"Reference base - pos {tmp_pos}, base {mut_change[1]} does not match the reference provided by user - base {cds_ref}")
                # track concecutive changes
                if len(concecutive_snp) == 0:
                    concecutive_snp.append(cds_pos)
                    combined_snp+=mut_change[2]
                else:
                    if cds_pos == concecutive_snp[-1] +1:
                        # update position in concec_pos
                        combined_snp+=mut_change[2]
                        # update basesin combined_bases
                        concecutive_snp.append(cds_pos)
                    elif cds_pos == concecutive_snp[-1] +2:
                        # get middle ref
                        m_ref = self._cds[concecutive_snp[-1]+1-1]
                        concecutive_snp.append(cds_pos)
                        combined_snp += m_ref+mut_change[2]
                    elif cds_pos > concecutive_snp[-1] +2:
                        # convert things in concecutive_snp and combined_snp into hgvs
                        hgvs = snp_to_hgvs(concecutive_snp, combined_snp, self._cds)
                        # update combined_snp and hgvs with current mut
                        concecutive_snp = [cds_pos]
                        combined_snp = mut_change[2]
                        mutations.append(hgvs)
        if len(concecutive_snp) !=0:
            hgvs = snp_to_hgvs(concecutive_snp, combined_snp, self._cds)
            mutations.append(hgvs)
        if len(delins) != 0:
            hgvs = delins_to_hgvs(self._cds, delins)
            mutations.append(hgvs)
        if len(mutations) == 1:
            mutations = f"c.{mutations[0]}"
        elif len(mutations) == 0:
            return [], []
        else:
            joined = ";".join(mutations)
            mutations = f"c.[{joined}]"
        return mutations, outside_mut

def snp_to_hgvs(concec_pos, combined_bases, cds):
    """
    helper function to obtain hgvs sstring given a list of positions and a string of bases combined.

    """
    if len(concec_pos) == 1:
        ref_base = cds[concec_pos[0]-1]
        hgvs = f"{concec_pos[0]}{ref_base}>{combined_bases}"
    else:
        hgvs = f"{concec_pos[0]}_{concec_pos[-1]}delins{combined_bases}"
    return hgvs


def delins_to_hgvs(cds_seq, delins):
    """
    convert a list of deletion and insertions (within 2 bp from eachother) to hgvs
    output hgvs should be ***_***delins***
    """
    # sort delins
    delins.sort()
    # convert delins to df
    print(delins)
    delins_df = pd.DataFrame(delins, columns=["pos", "base", "type", "qual"])
    types = delins_df["type"].unique()

    if len(types) == 1:
        # only one type of mutations in the list
        # means the list has len = 1
        delins = delins[0]
        if types[0] == "ins": # insertion
            hgvs = f"{delins[0]-1}_{delins[0]}ins{delins[1]}"
        else: # deletion
            # if it is a one bp deletion
            if len(delins[1]) == 1:
                hgvs = f"{delins[0]}del"
            # it is more than one bp
            else:
                hgvs = f"{delins[0]}_{delins[0]+len(delins[1])-1}del"

    else: # means that the len is >1
        start_pos = delins[0][0]
        end_pos = delins[-1][0]
        prev_pos = 0
        modified = ""
        for index, row in delins_df.iterrows():
            if index == 0:
                # if the first change is deletion
                if row["type"] == "del":
                    modified += cds_seq[row["pos"]:row["pos"]+len(row["base"])]
                    prev_pos += row["pos"] + len(row["base"])
                else:
                    modified += row["base"]
                    prev_pos += row["pos"]-1
            else:
                if row["type"] == "del":
                    # first we need to get whatever is between this change and previous change
                    mid = cds_seq[prev_pos: row["pos"]-len(row["base"])]
                    modified += mid
                    prev_pos += len(mid)
                else:
                    mid = cds_seq[prev_pos: row["pos"]]
                    modified += mid+row["base"]
                    prev_pos += len(mid)

        hgvs = f"{start_pos}_{end_pos}delins{modified}"
        #print("del+ins, ",hgvs)
    return hgvs

if __name__ == "__main__":

    # test delins_to_hgvs
    cds_seq = "CATCTT"
    print(cds_seq[1:3])
    delins = [[2, "A", "del"], [4, "GG", "ins"], [5, "T", "del"]]
    #delins = [[2, "G", "ins"], [4, "C", "del"]]
    delins_to_hgvs(cds_seq, delins)
