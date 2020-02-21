#~/lib/Python-3.6.4/python

## Read sam file (R1 and R2)
# count mutations in sam files 
# output mutation counts

import pandas as pd
import numpy as np
import re
import os
import sys
import glob
import logging
import argparse
import datetime
#from itertools import izip
from pathlib import Path

# modules in package 
import help_functions
import locate_mut
import locate_mut_tmp
import posterior

class readSam(object):

	def __init__(self, sam_r1, sam_r2, seq_lookup, tile_map, region_map, samples, output_dir, qual_filter, logf, log_level):
		"""
		sam_R1: read one of the sample
		sam_R2: read two of the sample
		seq_lookup: df contains DNA sequence, Protein sequence and positions mapped to each other

		tile_map: tile start and end pos
		region_map: region start and end pos

		samples: df with sample names and corresponding conditions, tiles
		output_dir: main output directory

		qual_filter: quality fileter number to sam file

		log level: settings for logging

		"""
		self._r1 = sam_r1
		self._r2 = sam_r2

		self._seq_lookup = seq_lookup

		self._tile_map = tile_map
		self._region_map = region_map

		self._qual = qual_filter

		self._output_counts_dir = output_dir

		# create logging file
		#self._log_dir = log_dir

		# sample information
		self._sample_id = os.path.basename(sam_r1).split("_")[0]
		self._sample_info = samples[samples["Sample ID"] == self._sample_id]

		self._sample_tile = self._sample_info["Tile ID"].values[0]
		self._tile_begins = (self._tile_map[self._tile_map["Tile Number"] == self._sample_tile]["Start AA"].values[0] *3)-2# beginning position of this tile (cds position)
		self._tile_ends = self._tile_map[self._tile_map["Tile Number"] == self._sample_tile]["End AA"].values[0] *3 # ending position of this tile (cds position)

		self._sample_condition = self._sample_info["Condition"].values[0]
		self._sample_tp = self._sample_info["Time point"].values[0]
		self._sample_rep = self._sample_info["Replicate"].values[0]

		if "_ds" in sam_r1:
			#log_f = "sample_"+ str(self._sample_id)+"_ds_mut_count.log"
			logging.basicConfig(filename=logf, filemode="a", format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p", level = log_level.upper())
			self._sample_counts_f = os.path.join(self._output_counts_dir,f"counts_sample{self._sample_id}_ds.csv")
		else:
			#log_f = "sample_"+ str(self._sample_id)+"_mut_count.log"
			logging.basicConfig(filename=logf, filemode="a", format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p", level = log_level.upper())

			self._sample_counts_f = os.path.join(self._output_counts_dir,f"counts_sample{self._sample_id}.csv")

		logging.info(f"Counting mutations in sample-{self._sample_id}")
		logging.info(f"Sam file input R1:{sam_r1}")
		logging.info(f"Sam file input R2:{sam_r2}")

		output_csv = open(self._sample_counts_f, "w")
		# write log information to counts output
		output_csv.write(f"#Sample:{self._sample_id}\n#Tile:{self._sample_tile}\n#Tile Starts:{self._tile_begins}\n#Tile Ends:{self._tile_ends}\n#Condition:{self._sample_condition}\n#Replicate:{self._sample_rep}\n#Timepoint:{self._sample_tp}\n""")
		output_csv.close()

	def _convert_sam(self, input_sam):
		"""
		trim sam file, convert it to a df with the following colmns:
		read_name, mapped_name, pos_start, mapQ, CIGAR, seq

		input sam files are filtered at this stage:
		* mapeed name should match input fasta sequence name
		* mapQ > mapQ cut off value
		* mdz indicates mutations in the read

		"""
		if "_R1_" in input_sam:
			read = "R1"
		else:
			read = "R2"
		read_counts = 0
		read_nomut = 0
		un_map = 0
		low_qual = 0
		trimmed = []
		# read sam file
		with open(input_sam, "r") as sam_file:
			for line in sam_file:
				if line.startswith("@"): # skip header line
					continue
				read_counts += 1
				line = line.split("\t")
				# create table to save information from sam file
				read_name = line[0]

				mapped_name = line[2]
				if mapped_name == "*": # if the read didn't map
					un_map +=1
					continue

				pos_start = line[3]

				mapQ = int(line[4])
				if mapQ < self._qual: # remove reads with mapQ < quality score filter
					low_qual +=1
					continue

				CIGAR = line[5]
				seq = line[9]
				quality = line[10]

				mdz = [i for i in line if "MD:Z:" in i]
				if len(mdz) != 0:
					mdz = mdz[0].split(":")[-1]
				else:
					mdz = ""
				if (not re.search('[a-zA-Z]', mdz)) and ("I" not in CIGAR):
					# if MDZ string only contains numbers 
					# and no insertions shown in CIGAR string
					# means there is no mutation in this read
					read_nomut +=1
					# remove reads that have no mutations in MDZ
					continue
				trimmed.append([read_name, mapped_name, pos_start, mapQ, CIGAR, mdz, seq, quality])

		trimmed = pd.DataFrame(trimmed, columns=["read_name", "mapped_name", "pos_start", "mapQ", "CIGAR", "mdz","seq", "qual"])
		output_csv = open(self._sample_counts_f, "a")
		output_csv.write(f"#Raw read depth for {read}:{read_counts}\n")
		output_csv.write(f"#Number of {read} reads without mutations:{read_nomut}\n#Number of {read} reads with qual < {self._qual}:{low_qual}\n#Number of {read} reads did not map to gene:{un_map}\n")
		output_csv.close()
		logging.info(f"Raw sequencing depth for {input_sam} (before filtering): {read_counts}")
		logging.info(f"Number of reads without mutations:{read_nomut}")
		return trimmed

	def _apply_thresh(self, mut_df):
		"""
		mut_df: dataframe contains mutation calls from R1 and R2
		return a list of mutations that passed user-defined threshold
		"""
		# group the mutations by position
		mut_group = mut_df.groupby("pos")
		for p in mut_group.groups.keys():
			n_mut = mut_group.get_group(p).shape[0]
			muts = mut_group.get_group(p)
			print(muts)
			if n_mut == 1:
				# in this case, only one read has the mutation
				basecall = muts.alt.tolist()
				phred = muts.qual.tolist()
				wt = muts.ref.tolist()[0]
				print(wt)
				post_prob = posterior.bayesian_variant_call(basecall, phred, wt, mut_rate=0.0025)
				print(post_prob)

	def _merged_main(self, full_seq, cds_seq):
		"""
		Read sam files at the same time, store mutations that passed filter
		"""
		read_pair = 0 # total pairs 
		un_map = 0 # total number of unmapped reads
		read_nomut = 0 # read pairs that have no mutations

		row = {}
		with open(self._r1, "r") as r1_f, open(self._r2, "r") as r2_f:
			for line_r1, line_r2 in zip(r1_f, r2_f): # this assumes both files have the same line #
				if line_r1.startswith("@") or line_r2.startswith("@"): # skip header lines
					continue
				#print(line_r1, line_r2)
				read_pair += 1

				line_r1 = line_r1.split("\t")
				line_r2 = line_r2.split("\t")

				mapped_name_r1 = line_r1[2]
				mapped_name_r2 = line_r2[2]
				if mapped_name_r1 == "*" or mapped_name_r2 == "*": # if one of the read didn't map to ref
					un_map +=1
					continue

				# check if read ID mapped 
				read_name_r1 = line_r1[0]
				read_name_r2 = line_r2[0]
				if read_name_r1 != read_name_r2: # reads are not in pairs
					logging.warning(f"Read ID did not map: {read_name_r1}, {read_name_r2}")

				# get starting position for r1 and r2
				pos_start_r1 = line_r1[3]
				pos_start_r2 = line_r2[3]

				# get CIGAR string
				CIGAR_r1 = line_r1[5]
				seq_r1 = line_r1[9]
				quality_r1 = line_r1[10]

				CIGAR_r2 = line_r2[5]
				seq_r2 = line_r2[9]
				quality_r2 = line_r2[10]

				mdz_r1 = [i for i in line_r1 if "MD:Z:" in i]
				mdz_r2 = [i for i in line_r2 if "MD:Z:" in i]

				if len(mdz_r1) != 0:
					mdz_r1 = mdz_r1[0].split(":")[-1]
				else:
					mdz_r1 = ""

				if len(mdz_r2) != 0:
					mdz_r2 = mdz_r2[0].split(":")[-1]
				else:
					mdz_r2 = ""

				if ((not re.search('[a-zA-Z]', mdz_r1)) and ("I" not in CIGAR_r1)) and ((not re.search('[a-zA-Z]', mdz_r2)) and ("I" not in CIGAR_r2)):
					# if MDZ string only contains numbers 
					# and no insertions shown in CIGAR string
					# means there is no mutation in this read
					# if both reads have no mutations in them, skip this pair
					read_nomut +=1
					# remove reads that have no mutations in MDZ
					continue

				# make the reads in the format of a dictionary
				# columns=["mapped_name", "pos_start", "qual", "CIGAR", "mdz","seq"])

				#row["mapped_name_r1"] = mapped_name_r1
				row["pos_start_r1"] = pos_start_r1
				row["qual_r1"] = quality_r1
				row["cigar_r1"] = CIGAR_r1
				row["mdz_r1"] = mdz_r1
				row["seq_r1"] = seq_r1

				#row["mapped_name_r2"] = mapped_name_r2
				row["pos_start_r2"] = pos_start_r2
				row["qual_r2"] = quality_r2
				row["cigar_r2"] = CIGAR_r2
				row["mdz_r2"] = mdz_r2
				row["seq_r2"] = seq_r2

				# pass this dictionary to locate mut
				# mut = locate_mut_main()
				# add mutation to mut list
				mut_parser = locate_mut_tmp.MutParser(row, full_seq, cds_seq, self._seq_lookup, self._tile_begins, self._tile_ends, logging)
				mut_parser._main()
				print(row)
				#break

		output_csv = open(self._sample_counts_f, "a")
		output_csv.write(f"#Raw read depth:{read_pair}\n")
		output_csv.write(f"#Number of read pairs without mutations:{read_nomut}\n#Number of reads did not map to gene:{un_map}\n")
		output_csv.close()
		logging.info(f"Raw sequencing depth: {read_pair}")
		logging.info(f"Number of reads without mutations:{read_nomut}")
		print(f"{read_nomut}, {read_pair}, {un_map}")


	def _main(self, full_seq, cds_seq):
		"""
		"""
		r1_df = self._convert_sam(self._r1)
		r2_df = self._convert_sam(self._r2)
		output_csv = open(self._sample_counts_f, "a")
		output_csv.write(f"#Aligned reads in R1:{r1_df.shape[0]}\n#Aligned reads in R2:{r2_df.shape[0]}\n")

		logging.info(f"Sequencing depth (after filtering) for R1: {r1_df.shape[0]}")
		logging.info(f"Sequencing depth (after filtering) for R2: {r2_df.shape[0]}")

		# merge r1 and r2 based on read ID
		merge = pd.merge(r1_df, r2_df, how="left", on="read_name", suffixes=("_r1", "_r2"))
		merge = merge.dropna()
		logging.info(f"After merging & filtering R1 sam file and R2 sam file, {merge.shape[0]} read-pairs remained for analysis")
		output_csv.write(f"#Final read-depth:{merge.shape[0]}\n")
		output_csv.write(f"#Comment: Final read-depth = Read pairs that passed the posterior threshold  with the same mutations\n")

		# facts['pop2050'] = facts.apply(lambda row: final_pop(row['population'],row['population_growth']),axis=1)
		mut_reads=0
		hgvs_output = []
		off_mut = []
		for index, row in merge.iterrows():
			# for each read pair in sam file, identify the mutation
			# process R1 and R2 together
			mut_parser = locate_mut.MutParser(row, full_seq, cds_seq, self._seq_lookup, self._tile_begins, self._tile_ends, logging)
			mp_update = mut_parser._get_seq()
			#mp_update = mut_parser._build_lookup()

			# get mutation from R1
			r1_mut =  mp_update._parse_mut(mp_update._r1_cigar, mp_update._r1_mdz, mp_update._r1_ref, mp_update._r1_read, mp_update._r1_pos, mp_update._r1_qual)
			#r1_mut["read"] = "R1"

			# get mutation from R2
			r2_mut = mp_update._parse_mut(mp_update._r2_cigar, mp_update._r2_mdz, mp_update._r2_ref, mp_update._r2_read, mp_update._r2_pos, mp_update._r2_qual)
			#r2_mut["read"]= "R2"
			# get overlap of mutations in R1 and mutations in R2
			#all_mut = pd.concat([r1_mut, r2_mut])
			# group the mutations by position
			#exit(0)
			both_mut = list(set(r1_mut) & set(r2_mut))
			off = ""
			if len(both_mut) !=0:
				both_mut.sort() # sort mutations based on the position
				# test hgvs
				hgvs, off= mp_update._get_hgvs(both_mut)
				if len(hgvs)!=0: # remove the case where mutations are not within the tile range
					hgvs_output.append(hgvs)
			if off != "":
				off_mut += off # mutations that are outside of the targeted tile

		off_mut = list(set(off_mut))
		output_csv.write(f"#Mutation positions outside of the tile: {off_mut}\n")
		output_csv.close()
		# convert list to df with one col
		hgvs_df = pd.DataFrame({"HGVS":hgvs_output})
		hgvs_counts = hgvs_df.HGVS.value_counts().reset_index()
		hgvs_counts.columns = ["HGVS", "count"]
		hgvs_counts.to_csv(self._sample_counts_f, mode="a", index=False)

		output_csv.close()


if __name__ == "__main__":

		parser = argparse.ArgumentParser(description='TileSeq mutation counts (for sam files)')
		parser.add_argument("-r1", "--read_1", help="sam file for R1", required=True)
		parser.add_argument("-r2", "--read_2", help="sam file for R2", required=True)
		parser.add_argument("-qual", "--quality", help="sam file mapQ filter", default=20)
		parser.add_argument("-o", "--output", help="Output folder", required = True)
		parser.add_argument("-logf", "--logf", help="Log file for this run", required=True)
		parser.add_argument("-log", "--log_level", help="Set log level: debug, info, warning, error, critical.", default = "debug")
		parser.add_argument("-p", "--param", help="csv paramter file", required = True)
		parser.add_argument("-t", "--thresh", help="posterior probability threshold")

		args = parser.parse_args()

		sam_r1 = args.read_1
		sam_r2 = args.read_2
		qual_filter = args.quality

		out = args.output
		param = args.param

		print("testing")
		# conver the csv file to json 
		# csv2json = os.path.abspath("src/csv2json.R")
		param_json = param.replace(".csv", ".json")
		#convert = f"Rscript {csv2json} {param} -o {param_json} -l stdout"
		#os.system(convert)

		# process the json file 
		project, seq, cds_seq, tile_map, region_map, samples = help_functions.parse_json(param_json)
		# build lookup table
		lookup_df = help_functions.build_lookup(seq.cds_start.item(), seq.cds_end.item(), cds_seq)

		# initialize the object
		MutCounts = readSam(sam_r1, sam_r2, lookup_df, tile_map, region_map, samples, out, qual_filter, args.logf, args.log_level)

		MutCounts._merged_main(seq, cds_seq)


