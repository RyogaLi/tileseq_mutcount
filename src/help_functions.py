#!~/lib/bin/bin/python3.6

# Helper functions
# 1. downsample fastq files into n reads
# 2. loading process animation function
# 3. parse input json file 

import pandas as pd
import random
import os
import sys
import json
import time
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Seq import MutableSeq

def downsample(n, r1, r2, output_path):
	"""
	Randomly downsample r1 and r2 to n reads
	"""

	record_number = 0

	base_r1 = os.path.basename(r1)
	base_r2 = os.path.basename(r2)

	output_r1 = output_path+"/"+base_r1.replace(".fastq", "_ds.fastq")
	output_r2 = output_path+"/"+base_r2.replace(".fastq", "_ds.fastq")

	with open(r1, "r") as i:
		r1_lines = sum([1 for line in i])
	with open(r2, "r") as j:
		r2_lines = sum([1 for line in j])

	total_records_r1 = int(r1_lines/4)
	total_records_r2 = int(r2_lines/4)

	if n > total_records_r1 or n > total_records_r2:
		os.system("cp "+r1+ " "+output_r1)
		os.system("cp "+r2+ " "+output_r2)
		return None
	else:
		if total_records_r1 <= total_records_r2:
			records_to_keep = set(random.sample(range(total_records_r1 + 1), n))
		else:
			records_to_keep = set(random.sample(range(total_records_r2 + 1), n))

	with open(r1, "r") as i_1:
		with open(output_r1, "w") as o1:

			for line1 in i_1:
				line2 = i_1.readline()
				line3 = i_1.readline()
				line4 = i_1.readline()
				if record_number in records_to_keep:
					o1.write(line1)
					o1.write(line2)
					o1.write(line3)
					o1.write(line4)
				record_number += 1


	record_number = 0
	with open(r2, "r") as i_2:
		with open(output_r2, "w") as o2:
			for line1 in i_2:
				line2 = i_2.readline()
				line3 = i_2.readline()
				line4 = i_2.readline()
				if record_number in records_to_keep:
					o2.write(line1)
					o2.write(line2)
					o2.write(line3)
					o2.write(line4)
				record_number += 1

	return output_r1, output_r2

def parse_csv(param_file, output_json, log):
	"""
	Convert user input cds file into json file
	Using	csv2json.R
	param_file
	output_json
	log
	"""
	pass


def parse_json(json_file):
	"""
	parse input json file
	return information about the run
	"""
	with open(json_file, "r") as jsonf:
		data = json.load(jsonf)
		project = data["project"] # project name
		seq = pd.DataFrame(data["template"], index=[0])
		cds_seq = seq.seq[0]
		cds_seq = cds_seq[int(seq.cds_start)-1: int(seq.cds_end)-1]
		# create dictionary to map tile/region # to start,end positions
		tile_map = pd.DataFrame(data["tiles"])
		region_map = pd.DataFrame(data["regions"])

		samples = pd.DataFrame(data["samples"])

	return project, seq, cds_seq, tile_map, region_map, samples


def loadingAnimation(process):
	"""
	loading animation to show when loading a process
	"""

	while process.is_alive():
		chars = "/—\|"
		for char in chars:
			sys.stdout.write('\r'+'downsamping reads ... '+char)
			time.sleep(.1)
			sys.stdout.flush()
	sys.stdout.write("\n")


def build_lookup(cds_start, cds_end, cds_seq):
    """
    build a lookup df for a given gene
    return a df with columns:
    temp_pos, na, dna_pos, protein, protein_pos
    """
    lookup_table = {}
    # list of template pos 
    temp = range(cds_start, cds_end)
    # list of coding DNA bases
    cDNA = [i for i in cds_seq]
    cDNA_pos = range(1, len(cDNA)+1)

    # list of protein bases
    protein = [i for i in Seq(cds_seq).translate()]
    protein_pos = range(1, len(protein)+1)
    lookup_table["temp_pos"] = temp
    lookup_table["cds"] = cDNA
    lookup_table["cds_pos"] = cDNA_pos
    lookup_table["protein"] = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in protein))
    lookup_table["protein_pos"] = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in protein_pos))

    df = pd.DataFrame(lookup_table)
    return df
