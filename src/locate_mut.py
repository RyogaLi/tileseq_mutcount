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

	def __init__(self, row, full_seq, cds_seq, seq_lookup, tile_s, tile_e, logging):
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

		self._tile_begins = tile_s
		self._tile_ends = tile_e
		# Note all the position in this df are 1-based
		self._seq_lookup = seq_lookup

		# get paired reads and their cigar scores
		self._reads = row

		# logger
		self._logging = logging

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
			#print(map_pos)

			if ins and map_pos >= ins[0]:
				inserted_pos += ins[1]
				ins = next(iter_ins, None)

			if "^" not in base:
				# this means a single nt change
				read_pos += match_len 
				mut_list.append(str(pos+read_pos-clip)+"|"+base+"|"+read[read_pos+inserted_pos-deleted_len])
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
				mut_list.append(str(pos+read_pos-clip)+"|"+base+"|"+read[read_pos-deleted_len])
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
		mut_list mutation list in the format
		translate a list of mut to aa changes
		return list of nt changes represent by hgvs strings
		"""
		# go through mutation list generated by _parse_mdz and _parse_cigar
		# make changes to the DNA sequence 
		# output protein changes for this read
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
				cds_pos = self._seq_lookup[self._seq_lookup.temp_pos == tmp_pos].cds_pos.item()
			except:
				continue

			if cds_pos < self._tile_begins or cds_pos > self._tile_ends:
				self._logging.warning(f"mutation at pos {cds_pos} which is not within the tile")
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
					if (mut_change[0] <= prev_pos+2) and (mut_change[2] != prev_change): # within 2bp of previous change
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
	# convert delins to df
	delins_df = pd.DataFrame(delins, columns=["pos", "base", "type"])
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
	return hgvs

if __name__ == "__main__":

	# test delins_to_hgvs
	cds_seq = "CATCTT"
	print(cds_seq[1:3])
	delins = [[2, "A", "del"], [4, "GG", "ins"], [5, "T", "del"]]
	#delins = [[2, "G", "ins"], [4, "C", "del"]]
	delins_to_hgvs(cds_seq, delins)
