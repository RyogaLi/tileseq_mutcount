#~/lib/Python-3.6.4/python

# Main script for sequencing analysis
# Input of this script: 
#		- JSON file containing input parameters
#   - Fastq files (gzipped)

# what does this script do?
#	
# 1. Read input csv file and convert to json (using csv2Json.R (you need to install required packages in R)
# 1. Read paramters from json file
# 2. Read input fastq file
#   2.a. Ramdonly select 30k reads and save a copy of downsampled fastq files 

# 3. Align fastq file to reference sequence, generate sam files 
# 4. From sam files, count mutations 
# 5. Output mutation counts to summary.csv

# modules
import pandas as pd
import numpy as np
import os
import glob
import argparse
import subprocess
import sys
import json
import logging
import datetime
import random
import sys, time, threading

# pakage modules
import settings
import alignment
import help_functions
import count_mut
import cluster


class MutCount(object):

	def __init__(self, param, project, seq, cds_seq, tile_map, region_map, samples, fastq_path, output_folder, log_level, env, n, skip=False):
		"""
		Initialize mutation counts
		Load input json file and parameters for this run
		project, seq, cds_seq, tile_map, region_map, samples
		output_folder: user input folder for storing time stamped output folders
		From param.json load all the parameters, save them into df
		"""
		self._param = param # parameter json file
		self._fastq_path = fastq_path
		self._fastq_list = glob.glob(fastq_path+"/*.fastq.gz")
		self._env = env
		self._n_reads = n

		# load parameter json file 
		self._project = project # project name

		# check if output folder exists
		if not os.path.isdir(output_folder):
			print(f"Output directory not found: {output_folder}")
			exit(1)

		if not skip:
			# get time stamp for this object
			time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

			# make time stamped output folder for this project
			self._output = self._project.replace(" ", "-")
			self._output = os.path.join(output_folder, self._output+ "_" +time)

			os.makedirs(self._output)
			logging.basicConfig(filename=os.path.join(self._output, "main.log"), filemode="w", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p", level = log_level.upper())
			logging.info(f"Analyzing {self._project} on {self._env}")
			logging.info(f"Fastq files are read from {self._fastq_path}")

		else:
			self._output = out
			# create main log file in this output folder
			# or append to the existing log file 
			logging.basicConfig(filename=os.path.join(self._output, "main.log"), filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p", level = log_level.upper())
			logging.info(f"Analyzing {self._project} on {self._env}")
			logging.info(f"Sam files are read from {self._output}")

		self._seq = seq
		self._cds_seq = cds_seq
		# ## validation: check if cds region is a multiple of 3
		# if len(self._cds_seq) % 3 != 0:
		#		logging.warning("CDS region lenth is not a multiple of 3")

		self._tile_map = tile_map
		self._region_map = region_map
		self._samples = samples

	def _align_sh_(self):
		"""
		1. downsample fastq files into n reads (specified by user, default 30,000)
		2. Align each fastq file to reference (submit one job for each alignment)
		3. Align each downsampled fastq file to reference (submit one job for each alignment)

		Output sam files into output folders (sam_files, ds_sam_files)

		"""
		## VALIDATE: if all samples are present in fastq files 
		# also check if any fastq file is empty

		sample_names = self._samples["Sample ID"].tolist()
		fastq_sample_id = [os.path.basename(i).split("_")[0] for i in self._fastq_list]
		fastq_sample_id = list(set(fastq_sample_id))

		## validation
		# check if input fastq files contains all the samples in paramter file
		if set(sample_names).issubset(fastq_sample_id):
			logging.info("All samples found")
		else:
			logging.error("fastq files do not match input samples.")
			logging.error("Program terminated due to error")
			exit()

		# create mapping for R1 and R2 for each sample
		# ONLY samples provided in the parameter files will be analyzed
		fastq_map = []
		for r1 in self._fastq_list:
			ID = os.path.basename(r1).split("_")[0]
			if ("_R1_" in r1) and (ID in sample_names):
				r2 = r1.replace("_R1_", "_R2_")
				fastq_map.append([r1,r2])

		# convert to df
		fastq_map = pd.DataFrame(fastq_map, columns=["R1", "R2"])

		### make folder to store all downsampled fastq files
		ds_output = os.path.join(self._output, "ds_files")
		os.system("mkdir "+ds_output)
		## downsample all the files into dir
		logging.info(f"Downsampling fastq files to {self._n_reads} reads. ")
		loading_process = threading.Thread(name = "process", target=ds_process, args = (fastq_map, self._n_reads, ds_output))
		loading_process.start()
		help_functions.loadingAnimation(loading_process)
		loading_process.join()
		logging.info("Downsampling finished!")

		# make folder to store alignment sam files 
		sam_output = os.path.join(self._output, "sam_files/")
		os.system("mkdir "+sam_output)

		# make folder to sotre downsampled alignment sam files
		ds_sam_output = os.path.join(self._output, "ds_sam_files/")
		os.system("mkdir "+ds_sam_output)

		# make folder to store ref
		ref_path = os.path.join(self._output, "ref/")
		os.system("mkdir "+ ref_path)

		# GURU
		if self._env == "GURU":
			# make sh files to submit to GURU
			sh_output = os.path.join(self._output, "GURU_sh")
			os.system("mkdir "+sh_output)

			# for each pair of fastq files, submit jobs to run alignment 
			# submit alignments for both full run and downsampled run 
			# for the alignment module
			# it takes the following arguments: R1, R2, ref, Output sam
			# the output log (bowtie log) would be in the same dir
			logging.info("Writing sh files for alignment (GURU)")
			cluster.alignment_sh_guru(fastq_map, self._project, self._seq.seq.item(), ref_path, sam_output, ds_sam_output, sh_output, logging)
			logging.info("Alignment jobs are submitted to GURU..")

			# wait for alignment to finish and call mutations

		elif self._env == "BC2":
			# get current user name 
			cmd = ["whoami"]
			process = subprocess.run(cmd, stdout=subprocess.PIPE)
			userID = process.stdout.decode("utf-8").strip()

			# make sh files to submit to BC
			sh_output = os.path.join(self._output, "BC_sh")
			os.system("mkdir "+sh_output)

			logging.info("Writing sh files for alignment (BC2)")
			sam_df = cluster.alignment_sh_bc2(fastq_map, self._project, self._seq.seq.item(), ref_path, sam_output, ds_sam_output, sh_output, logging)
			logging.info("Alignment jobs are submitte to BC2. Check pbs-output for STDOUT/STDERR")

			# get number of jobs running
			check = ["qstat", "-u", userID]
			check_process = subprocess.run(check, stdout=subprocess.PIPE)
			n_jobs = check_process.stdout.decode("utf-8").strip()
			while n_jobs != '':
				n = n_jobs.count("\n")
				logging.info(f"{n-2} jobs running ....")
				# wait for 10 mins
				time.sleep(300)
				check_process = subprocess.run(check, stdout=subprocess.PIPE)
				n_jobs = check_process.stdout.decode("utf-8").strip()

			return sam_df

	def _mut_count(self, sam_df=pd.DataFrame()):
		"""
		Count mutations in input sam files
		1. Log path to sam files and downsampled sam files
		2. If sam df not provided, make df for all the sam files
		3. For each pair of sam files, submit jobs to the cluster
		4. Log to main log
		"""
		# make dir for mut counts
		time_stamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
		mut_output_dir = os.path.join(self._output, time_stamp + "_mut_call")
		os.makedirs(mut_output_dir)

		log_dir = os.path.join(mut_output_dir, "mut_log")
		os.makedirs(log_dir)

		#log parameters for this run
		param_logf = os.path.join(mut_output_dir, "param.log")
		with open(param_logf, "w") as p:
			p.write("Input parameters for this run ...")
			p.write(f"Project name:")
			p.write(f"")


		if sam_df.empty: # skip alignment 
			logging.info(f"Skipping alignment...")
			logging.info(f"Analyzing sam files in {self._output}")
			# make sam_df 
			# get list of sam files in the ouput folder given by user 
			# folder that stores alignment sam files 
			sam_output = os.path.join(self._output, "sam_files/")
			sam_list = os.listdir(sam_output)

			# folder that stores downsampled alignment sam files
			ds_sam_output = os.path.join(self._output, "ds_sam_files/")
			ds_sam_list = os.listdir(ds_sam_output)

			# merge sam files and ds sam files into a df
			sam_df = []
			for f in sam_list:
				sample = os.path.basename(f).split("_")[0]
				if "_R1_" in f: # take the read one files 
					sam_r1 = os.path.join(sam_output, f)
					sam_r2 = os.path.join(sam_output, f.replace("_R1_", "_R2_"))
					# find the file in ds_sam_files
					sam_ds_r1 = glob.glob(f"{ds_sam_output}/{sample}*_R1_*.sam")[0]
					sam_ds_r2 = glob.glob(f"{ds_sam_output}/{sample}*_R2_*.sam")[0]

					sam_df.append([sam_r1, sam_r2, sam_ds_r1, sam_ds_r2])
			# convert sam_df to dataframe 
			# col names = r1_sam, r2_sam, r1_sam_ds, r2_sam_ds
			sam_df = pd.DataFrame.from_records(sam_df, columns = ["r1_sam", "r2_sam", "r1_sam_ds", "r2_sam_ds"])

		if self._env == "BC2":
			sh_output = os.path.join(self._output, "BC_sh")

		logging.info("Submitting mutation counts jobs to BC...")

		cluster.mut_count_sh_bc(sam_df, mut_output_dir, self._param, sh_output, log_dir, logging)

		# get number of jobs running
		cmd = ["whoami"]
		process = subprocess.run(cmd, stdout=subprocess.PIPE)
		userID = process.stdout.decode("utf-8").strip()

		check = ["qstat", "-u", userID]
		check_process = subprocess.run(check, stdout=subprocess.PIPE)
		n_jobs = check_process.stdout.decode("utf-8").strip()
		while n_jobs != '':
			n = n_jobs.count("\n")
			logging.info(f"{n-2} jobs running ....")
			# wait for 5 mins
			time.sleep(300)
			check_process = subprocess.run(check, stdout=subprocess.PIPE)
			n_jobs = check_process.stdout.decode("utf-8").strip()

		# job complete if nothing is running 
		# go through mutation call files generated and log files without any mutations
		logging.info("Check mutation counts file ...")
		mutcount_list = glob.glob(os.path.join(mut_output_dir, "counts_sample*.csv"))
		logging.info(f"{len(mutcount_list)} mutation counts file generated")
		for f in mutcount_list:
			mut_n = 0
			# double check if each output file has mutations 
			with open(f, "r") as mut_output:
				for line in mut_output:
					# skip header
					if "c." in line:
						mut_n += 1
			if mut_n == 0:
				logging.error(f"{f} has 0 variants! Check mut log for this sample.")
			else:
				logging.info(f"{f} has {mut_n} variants")
		logging.info("Job complete")


def ds_process(fastq_map, n, ds_output):
	"""
	fastq_map: dataframe with two columns, path to [R1, R2]
	return: dataframe with four columns, path to [R1, R2, R1_ds, R2_ds]
	Call downsample() in ds.py
	"""
	ds_r1files, ds_r2files = [], []
	for index, row in fastq_map.iterrows():
		os.system("gunzip "+row["R1"])
		os.system("gunzip "+row["R2"])

		unziped_r1 = row["R1"].replace(".gz", "")
		unziped_r2 = row["R2"].replace(".gz", "")

		r1_ds, r2_ds = help_functions.downsample(n, unziped_r1, unziped_r2, ds_output)
		os.system(f"gzip {unziped_r1}")
		os.system(f"gzip {unziped_r2}")
		r1_ds = r1_ds + ".gz"
		r2_ds = r2_ds + ".gz"
		ds_r1files.append(r1_ds)
		ds_r2files.append(r2_ds)
		#print(row["R1"])
		#break
	# gzip everything 
	os.system(f"gzip {ds_output}/*.fastq")
	fastq_map[["r1_ds", "r2_ds"]] = pd.DataFrame(list(zip(ds_r1files, ds_r2files)))
	return fastq_map


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='TileSeq mutation counts')
	parser.add_argument("-f", "--fastq", help="Path to all fastq files you want to analyze", required= True)
	parser.add_argument("-o", "--output", help="Output folder", required = True)
	parser.add_argument("-log", "--log_level", help="set log level: debug, info, warning, error, critical.", default = "debug")
	parser.add_argument("-p", "--param", help="csv paramter file", required = True)
	parser.add_argument("-env", "--environment", help= "The cluster used to run this script", default="BC2")
	parser.add_argument("-n", "--n_reads", help="Used for downsampling the files. n_reads will remain", default = 30000)
	parser.add_argument("--skip_alignment", action="store_true", help="skip alignment for this analysis, ONLY submit jobs for counting mutations in existing output folder")
	parser.add_argument("--skip_ds", action="store_true", help="skip downsampling of the fastq files")
	parser.add_argument("-thresh", "--threshold", help="Posterior threshold for filtering mutations")

	args = parser.parse_args()

	f = args.fastq
	out = args.output
	param = args.param
	env = args.environment
	n = args.n_reads

	log_level = args.log_level

	print(f"Log level: {log_level}")
	print(f"Convert {param} to json format")
	print("Output json file will be saved in the same dir as input parameter.csv file")

	cds2json = os.path.abspath("csv2json.R")
	param_json = param.replace(".csv", ".json")

	convert = f"Rscript {cds2json} {param} -o {param_json} -l stdout"
	os.system(convert)

	# read json file 
	project, seq, cds_seq, tile_map, region_map, samples = help_functions.parse_json(param_json)

	# check flag
	# if --skip-alignment, only submit jobs for mutation counts
	if args.skip_alignment:
		# Initialize MutCount main 
		mc = MutCount(param, project, seq, cds_seq, tile_map, region_map, samples, f, out, log_level, env, n, skip=True)
		print("skipping alignment ..")
		print(f"Analyzing sam files in {out}")
		mc._mut_count()
	else:
		# alignment
		# return job ID list (for checking if the jobs are running still)
		mc = MutCount(param, project, seq, cds_seq, tile_map, region_map, samples, f, out, log_level, env, n, skip=False)
		sam_df = mc._align_sh_()
		mc._mut_count(sam_df)
