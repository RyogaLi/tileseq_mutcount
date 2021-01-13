#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import os
import pandas as pd
import json
import argparse
from TileSeqMut.main import log

#for two runs merge the mutation counts files into output folder


def merge_runs(input_dir1, input_dir2, map_df, output):
	"""
	:param input_dir1: Directory 1 contains mutation counts files
	:param input_dir2: Directory 2 contains mutation counts files
	:param map_df: A dataframe that maps the samples to be merged
	:param output: output directory to save merged mut_counts
	"""
	# make a new file

	new_samples = []
	for index, row in map_df.iterrows():
		merged_header = ""
		mut_count_f1 = os.path.join(input_dir1, f"counts_sample_{row['Sample ID_run1']}.csv")
		mut_count_f2 = os.path.join(input_dir2, f"counts_sample_{row['Sample ID_run2']}.csv")
		new_samples.append(f"{row['Sample ID_run1']}-{row['Sample ID_run2']}")
		# join headers
		with open(mut_count_f1) as file1, open(mut_count_f2) as file2:
			for line1, line2 in zip(file1, file2):
				if "#" not in line1: break
				if "Final read-depth:" in line1:
					n_depth_1 = int(line1.split(":")[1])
					n_depth_2 = int(line2.split(":")[1])
					total_depth = n_depth_1+n_depth_2
					merged_header+= f"#Final read-depth:{total_depth}\n"
					continue
				merged_header += line1.strip() + "; " + line2
		output_f = open(os.path.join(output, f"counts_sample_{row['Sample ID_run1']}-{row['Sample ID_run2']}.csv"), "w")
		output_f.write(merged_header)
		# skip header and read the mutations
		df1 = pd.read_csv(mut_count_f1, skiprows=18)
		df2 = pd.read_csv(mut_count_f2, skiprows=18)
		merge = [df1, df2]
		merge = pd.concat(merge)
		merge_counts = merge.groupby(by=["HGVS"], as_index=False)["count"].sum().sort_values("count", ascending=False)
		# save header and counts to file
		merge_counts.to_csv(output_f, mode="a", index=False)
		output_f.close()

	return new_samples


def parse_params(param1, param2, output, main_log):
	"""
	Given two parameter sheets, find samples in common and output
	df for later merging
	:param param1: parameter sheet for the first run
	:param param2: parameter sheet for the second run
	:param main_log: logger
	:return map_df: df with two cols [samples_1, samples_2]
	"""
	update_json = {}
	with open(param1, "r") as p1:
		with open(param2, "r") as p2:
			param_one = json.load(p1)
			param_two = json.load(p2)
			samples1 = param_one["samples"]
			samples2 = param_two["samples"]
			samples1_df = pd.DataFrame(samples1)
			samples2_df = pd.DataFrame(samples2)
	# from these two files find samples in the same tile, with the same condition names
	# same replicates and tp
	# log how many samples will be merged
	merge = pd.merge(samples1_df, samples2_df, how="inner", on=["Tile ID", "Condition", "Time point", "Replicate"],
					 suffixes=["_run1", "_run2"])
	unique_Tid = merge["Tile ID"].drop_duplicates().tolist()
	unique_Conditions = merge["Condition"].drop_duplicates().tolist()
	merge["new_samples"] = merge["Sample ID_run1"] + "-" + merge["Sample ID_run2"]
	# save new sample names to output
	merge[["new_samples", "Tile ID", "Condition", "Time point", "Replicate"]].to_csv(os.path.join(output, "new_sample_names.csv"), index=False)
	main_log.info(f"{merge.shape[0]} samples will be merged")
	main_log.info(f"Tiles: {unique_Tid}")
	main_log.info(f"Condition: {unique_Conditions}")

	return merge


def merge_main(args):

	# check args and make corresponding output dir
	# log input/output information
	if not os.path.isdir(args.output):
		os.mkdir(args.output)
	logging = log(args.output, args.log_level)
	main_log = logging.getLogger("mergeRuns.log")
	main_log.info(f"Input dir 1: {args.dir1} with Parameter Sheet: {args.paramOne}")
	main_log.info(f"Input dir 2: {args.dir2} with Parameter Sheet: {args.paramTwo}")
	main_log.info(f"Merges mutation count file will be saved to {args.output}")

	map_df = parse_params(args.paramOne, args.paramTwo, args.output, main_log)
	new_samples = merge_runs(args.dir1, args.dir2, map_df, args.output)

	main_log.info(f"Merged files are saved to {args.output}")

	# make a json file with merged sample names (copy the other information from old json files)
	new_paramsheet = ""

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='TileSeq mutation counts')
	# user input arguments
	parser.add_argument("-p1", "--paramOne", help="Path to parameter sheet (JSON) for the first run", type=str, required=True)
	parser.add_argument("-p2", "--paramTwo", help="Path to parameter sheet (JSON) for the second run", type=str, required=True)
	parser.add_argument("-d1", "--dir1", help="Path to *_mut_count folder for the first run", type=str, required=True)
	parser.add_argument("-d2", "--dir2", help="Path to *_mut_count folder for the second run", type=str, required=True)
	parser.add_argument("-o", "--output", help="Output folder path (this folder will be made if not exists when "
											   "running the script",
						type=str, required=True)

	parser.add_argument("-log", "--log_level", help="set log level: debug, \
	        info, warning, error, critical. (default = debug)", type=str, default="debug")

	args = parser.parse_args()

	merge_main(args)
