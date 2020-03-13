#~/lib/Python-3.6.4/python

# Main script for sequencing analysis
# Input of this script: 
#		- JSON file containing input parameters
#	- Fastq files (gzipped)

# what does this script do?
#	
# 1. Read input csv file and convert to json (using csv2Json.R (you need to install required packages in R)
# 1. Read paramters from json file
# 2. Read input fastq file
#	2.a. Ramdonly select 30k reads and save a copy of downsampled fastq files 

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
import shutil
import sys, time, threading

# pakage modules
import settings
import alignment
import help_functions
import count_mut
import cluster


class MutCount(object):

	def __init__(self, param, fastq_path, output_folder, log_level, min_cover, mt, at, env, qual, logging):
		"""
		Initialize mutation counts
		Load input json file and parameters for this run
		project, seq, cds_seq, tile_map, region_map, samples
		output_folder: user input folder for storing time stamped output folders
		From param.json load all the parameters, save them into df
		param: JSON file
		fastq_path: path to fastq files
		output_folder: Path to output dir
		log_level: log level for the logging object
		min_cover: min cover percentage
		mt: mutation call time
		at: alignement time
		env: which cluster
		qual: posterior quality cutoff
		main_log: logging object
		"""
		self._param = param # parameter json file
		self._fastq_path = fastq_path
		self._fastq_list = glob.glob(fastq_path+"/*.fastq.gz")
		self._env = env
		self._cutoff = qual
		self._min_cover = min_cover
		self._mt = mt
		self._at = at
		self._output = output_folder
		# parse parameter json file
		self._project, self._seq, self._cds_seq, self._tile_map, self._region_map, self._samples = help_functions.parse_json(param)
		self._project = project.replace(" ", "_") # project name
		self._log = logging.getLogger("main.log")

	def _init_skip(self, skip, r1=None, r2=None):
		"""
		if skip == T: skip alignment
			make suboutput dir with time stamp
		else: not skip alignment
			make main output dir with time stamp
		"""
		if skip: # only run mutation counts on the aligned samples
				if r1 and r2:
					self._r1 = r1
					self._r2 = r2
					self._log.info("A pair of SAM files are provided:")
					self._log.info(f"SAM R1: {r1}")
					self._log.info(f"SAM R2: {r2}")
				else:
					self._log.info(f"Sam files are read from {self._output}")

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
				self._log.info("All samples found")
				self._log.info(f"In total there are {len(list(set(sample_names)))} samples in the csv file")
				self._log.info(f"In total there are {len(fastq_sample_id)} fastq files")
		else:
				self._log.error("fastq files do not match input samples.")
				self._log.error("Program terminated due to error")
				exit(1)

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

		# make folder to store alignment sam files 
		sam_output = os.path.join(self._output, "sam_files/")
		os.system("mkdir "+sam_output)

		# make folder to store ref
		ref_path = os.path.join(self._output, "ref/")
		os.system("mkdir "+ ref_path)

		# GURU
		if self._env == "GURU":
				# make sh files to submit to GURU
				sh_output = os.path.join(self._output, "GURU_sh")
				os.system("mkdir "+sh_output)

				# for each pair of fastq files, submit jobs to run alignment 
				# it takes the following arguments: R1, R2, ref, Output sam
				# the output log (bowtie log) would be in the same dir
				logging.info("Writing sh files for alignment (GURU)")
				cluster.alignment_sh_guru(fastq_map, self._project, self._seq.seq.item(), ref_path, sam_output, sh_output, logging)
				logging.info("Alignment jobs are submitted to GURU..")

				# wait for alignment to finish and call mutations

		elif self._env == "BC2" or self._env == "DC" or self._env == "BC":
				# get current user name 
				cmd = ["whoami"]
				process = subprocess.run(cmd, stdout=subprocess.PIPE)
				userID = process.stdout.decode("utf-8").strip()

				# make sh files to submit to BC
				sh_output = os.path.join(self._output, "BC_aln_sh")
				os.system("mkdir "+sh_output)

				if self._env == "BC2" or self._env == "BC":
						# make sh files to submit to BC
						sh_output = os.path.join(self._output, "BC_aln_sh")
						os.system("mkdir "+sh_output)

						self._log.info("Submitting alignment jobs to BC/BC2...")
						sam_df = cluster.alignment_sh_bc2(fastq_map, self._project, self._seq.seq.values.item(), ref_path, sam_output, sh_output, self._at, logging)
						self._log.info("Alignment jobs are submitte to BC2. Check pbs-output for STDOUT/STDERR")
				else:
						# make sh files to submit to DC
						sh_output = os.path.join(self._output, "DC_aln_sh")
						os.system("mkdir "+sh_output)

						self._log.info("Submitting alignment jobs to DC...")
						sam_df = cluster.alignment_sh_dc(fastq_map, self._project, self._seq.seq.values.item(), ref_path, sam_output, sh_output, self._at, logging)
						self._log.info("Alignment jobs are submitte to DC. Check pbs-output for STDOUT/STDERR")

				# get number of jobs running
				check = ["qstat", "-u", userID]
				check_process = subprocess.run(check, stdout=subprocess.PIPE)
				n_jobs = check_process.stdout.decode("utf-8").strip()
				n = n_jobs.count("\n")
				logging.info(f"{n-2} jobs running ....")
				while n_jobs != '':
						n = n_jobs.count("\n")
						logging.info(f"{n-2} jobs running ....")
						# wait for 10 mins
						time.sleep(300)
						check_process = subprocess.run(check, stdout=subprocess.PIPE)
						n_jobs = check_process.stdout.decode("utf-8").strip()
				logging.info(f"All alignment finished")
				# check how many sam files generated in the sam_files
				n_sam = len(os.listdir(sam_output))
				logging.info(f"{n_sam} sam files generated in {sam_output}")

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

		if sam_df.empty: # skip alignment 
				logging.info(f"Skipping alignment...")
				logging.info(f"Analyzing sam files in {self._output}")
				# make sam_df 
				# sam df is built by reading the input csv file
				sample_names = self._samples["Sample ID"].tolist()
				sam_dir = os.path.join(self._output, "sam_files/")
				sam_list = []
				for i in sample_names:
						sam_f = glob.glob(f"{sam_dir}{i}_*.sam") # assume all the sam files have the same name format (id_*.sam)
						if len(sam_f) < 2:
								logging.error(f"SAM file for sample {i} not found")
								exit(1)
						elif len(sam_f) == 2:
								for sam in sam_f:
										if "_R1_" in sam: # for read one
												sam_r1 = sam
										else:
												sam_r2 = sam
						sam_list.append([sam_r1, sam_r2])

				# convert sam_df to dataframe 
				sam_df = pd.DataFrame.from_records(sam_list, columns = ["r1_sam", "r2_sam"])

		if self._env == "BC2" or self._env == "DC" or self._env == "BC":
				sh_output = os.path.join(mut_output_dir, "BC_mut_sh")
				os.mkdir(sh_output)
				if self._env == "BC2" or self._env == "BC":
						logging.info("Submitting mutation counts jobs to BC2...")

						cluster.mut_count_sh_bc(sam_df, mut_output_dir, self._param, sh_output, self._min_cover, self._mt, log_dir, logging, self._cutoff)
						logging.info("All jobs submitted")
				else:
						logging.info("Submitting mutation counts jobs to DC...")

						cluster.mut_count_sh_dc(sam_df, mut_output_dir, self._param, sh_output, self._min_cover, self._mt, log_dir, logging, self._cutoff)
						logging.info("All jobs submitted")


				# get number of jobs running
				cmd = ["whoami"]
				process = subprocess.run(cmd, stdout=subprocess.PIPE)
				userID = process.stdout.decode("utf-8").strip()

				check = ["qstat", "-u", userID]
				check_process = subprocess.run(check, stdout=subprocess.PIPE)
				n_jobs = check_process.stdout.decode("utf-8").strip()
				n = n_jobs.count("\n")
				logging.info(f"{n-2} jobs running ....")
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
				all_tmp = os.path.join(mut_output_dir, "*_tmp.csv")
				merged = os.path.join(mut_output_dir, "info.csv")
				cmd = f"cat {all_tmp} > {merged}"
				os.system(cmd)
				#os.system(f"rm {all_tmp}")
				logging.info("Job complete")

	def _main(self):
		"""
		"""
		# submit jobs for alignment if self._align is true
		# obtain a list of running job id and check if the job finished 
		# if all the jobs submitted are finished, proceed to mutation counting
		# otherwise wait for 30min
		# return is sam_df


		# submit jobs for mutation counting 
		# 

def main(f, out, param, env, qual, min_cover, log_level, r1sam, r2sam, mt, at, run_name):
	"""
	Main for fastq2counts

	"""
	# get time stamp for this object
	time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

	# check if output dir exists
	if not os.path.isdir(out):
			print(f"Output directory not found: {output_folder}")
			exit(1)

	if not os.path.isdir(f):
			print(f"Fastq file path not found: {fastq_path}")
			exit(1)

	# try convert csv to json in the same dir as the csv file
	# convert csv file to json
	#param_path = shutil.copy(param, output_dir, *, follow_symlinks=True)
	if param.endswith(".csv"):
			param_json = param.replace(".csv", ".json")
			convert = f"Rscript {settings.CSV2JSON} {param} -o {param_json} --srOverride"
			os.system(convert)
			# check if json file exist
			if not os.path.isfile(param_json):
					print("Json file does not exist, check conversion!")
					exit(1)
	elif param.endswith(".json"):
			param_json = param
	else:
			print("Please provide valid paramter file format (csv or json)")
			exit(1)

	if not os.path.isfile(param_json):
			print("Json file does not exist, check conversion!")
			exit(1)
	# get basename for the param file
	param_base = os.path.basename(param_json)
	output_dir = ""
	if skip_alignment:
			# if we want to skip the alignment part
			# CASE 1: user provided a json file and R1sam + R2sam
			if r1sam and r2sam:
					# analyze one pair of r1 and r2 file
					if not r1sam.endswith(".sam") or not r2sam.endswith(".sam"):
							print("Please provide SAM files")
							exit(1)
					# load json file 
					param_path = os.path.join(out, param_base)
					if not os.path.isfile(param_path):
							param_path = shutil.copy(param, output_dir, *, follow_symlinks=True)
					#project, seq, cds_seq, tile_map, region_map, samples = help_functions.parse_json(param)
					# set up logging
					main_log = log(output_dir, log_level.upper(), "main")
					mc = mutcount(param_path, f, out, log_level, min_cover, mt, at, env, qual, main_log)
					mc._init_skip(skip=True, r1=r1sam, r2=r2sam)

			else:
					# user provided a csv file and path to sam files
					# for each pair of sam file we would submit one job to the cluster for mutation counting
					# make time stamped output folder for this project
					output_dir = os.path.join(out, run_name+ "_" +time+"_mut_count")
					os.makedirs(output_dir)
					# load json file 
					param_path = os.path.join(output_dir, param_base)
					if not os.path.isfile(param_path):
							param_path = shutil.copy(param, output_dir, *, follow_symlinks=True)
					mc = mutcount(param_json, f, output_dir, log_level, min_cover, mt, at, env, qual, main_log)
					mc._init_skip(skip=True)
					main_log = log(output_dir, log_level.upper(), "main")

	else:
			# alignment
			# make time stamped output folder for this project
			output_dir = os.path.join(output_folder, run_name+ "_" +time)
			os.makedirs(output_dir) # make directory to save this run 
			param_path = os.path.join(out, param_base)
			if not os.path.isfile(param_path):
					param_path = shutil.copy(param, output_dir, *, follow_symlinks=True)
			main_log = log(output_dir, log_level.upper(), "main")
			mc = mutcount(param_json, f, output_dir, log_level, min_cover, mt, at, env, qual, main_log)
			mc._init_skip(skip=False)


	if output_dir != "":
			# this parameter file will be saved in the main output dir 
			param_f = open(os.path.join(output_dir, "param.log"), "a")
			# log - time for this run (same as the output folder name)
			param_f.write(f"Run started: {time}\n")
			param_f.write(f"Run name: {run_name}\n")
			param_f.write(f"This run was on the cluster: {enc}\n")
			param_f.write(f"Input parameter file used: {param}\n")
			param_f.write(f"Posterior probability cutoff for mutation counts: {qual}\n")
			param_f.write(f"Min_ cover percentage: {min_cover}\n")
			param_f.close()

	# start the run 
	mc._main()

def log(output_dir, log_level.upper(), log_name):
	"""
	Make a logging object which writes to the main.log in output_dir
	"""
	# init main log
	logging.basicConfig(filename=os.path.join(output_dir, "main.log"),
					filemode="a",
					format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
					datefmt="%m/%d/%Y %I:%M:%S %p",
					level = log_level)
	# define a Handler which writes INFO messages or higher to the sys.stderr
	console = logging.StreamHandler()
	console.setLevel(logging.log_level)
	# set a format which is simpler for console use
	formatter = logging.Formatter('%(name)-8s: %(levelname)-4s %(message)s')
	# tell the handler to use this format
	console.setFormatter(formatter)
	# add the handler to the root logger
	logging.getLogger('').addHandler(console)

	return logging

if __name__ == "__main__":

		parser = argparse.ArgumentParser(description='TileSeq mutation counts')
		parser.add_argument("-f", "--fastq", help="Path to all fastq files you want to analyze", required= True)
		parser.add_argument("-o", "--output", help="Output folder", required = True)
		parser.add_argument("-log", "--log_level", help="set log level: debug, info, warning, error, critical. (default = debug)", default = "debug")
		parser.add_argument("-p", "--param", help="csv paramter file", required = True)
		parser.add_argument("-env", "--environment", help= "The cluster used to run this script (default = DC)", default="DC")
		parser.add_argument("--skip_alignment", action="store_true", help="skip alignment for this analysis, ONLY submit jobs for counting mutations in existing output folder")
		parser.add_argument("-qual", "--quality", help="Posterior threshold for filtering mutations (default = 0.99)", default = 0.99)
		parser.add_argument("-name", help="Name for this run", required = True)
		parser.add_argument("-min", "--min_cover", help="Minimal % required to cover the tile (default = 0.4)", default = 0.6)
		parser.add_argument("-at", type = int, help="Alignment time (default = 8h)", default = 8)
		parser.add_argument("-mt", type = int, help="Mutation call time (default = 36h)", default = 36)

		parser.add_argument("-r1", help="r1 SAM file")
		parser.add_argument("-r2", help="r2 SAM file")

		args = parser.parse_args()

		# convert all the args to variables
		f = args.fastq # path to fastq files
		out = args.output # path to ouptut dir
		param = args.param # path to parameter.csv
		env = args.environment # environment settings (DC, BC, BC2)
		qual = float(args.quality) # posterior quality cutoff
		min_cover = float(args.min_cover) # min coverage
		log_level = args.log_level # log level defaultl
		r1sam = args.r1 # sam file r1
		r2sam = args.r2 # sam file r2
		mt = args.mt # mutation count time required 
		at = args.at # alignment time required
		name = args.name
		print(" **** NOTE: Before you run this pipeline, please check settings.py to update program paths **** ")
		main(f, out, param, env, qual, min_cover, log_level, r1sam, r2sam, mt, at, name)
