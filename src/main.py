#~/.pyenv/shims/python3

# Main script for sequencing analysis
# Input of this script: 
#		- JSON file containing input parameters
#   - Fastq files (gzipped)

# what does this script do?
#
# 1. Read input JSON file, read paramters 
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
import glob
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


class MutCount(object):

	def __init__(self, project, seq, cds_seq, tile_map, region_map, samples, fastq_path, output_folder, log_level, env, n):
		"""
		Initialize mutation counts
		Load input json file and parameters for this run
		project, seq, cds_seq, tile_map, region_map, samples
		output_folder: user input folder for storing time stamped output folders
		From param.json load all the parameters, save them into df
		"""
		self._fastq_path = fastq_path
		self._fastq_list = glob.glob(fastq_path+"/*.fastq.gz")
		self._env = env
		self._n_reads = n

		# check if output folder exists
		if not os.path.isdir(output_folder):
			print(f"Output directory not found: {output_folder}")
			exit(1)

		# get time stamp for this object
		time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

		# load parameter json file 
		self._project = project # project name

		# make time stamped output folder for this project
		self._output = self._project.replace(" ", "-")
		self._output = os.path.join(output_folder, self._output+ "_" +time)

		os.makedirs(self._output)

		# create main log file in this output folder
		logging.basicConfig(filename=os.path.join(self._output, "main.log"), filemode="w", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p", level = log_level.upper())
		logging.info(f"Analyzing {self._project} on {self._env}")
		logging.info(f"Fastq files are read from {self._fastq_path}")

		self._seq = seq
		self._cds_seq = cds_seq
		## validation: check if cds region is a multiple of 3
		if len(self._cds_seq) % 3 != 0:
			logging.warning("CDS region lenth is not a multiple of 3")

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
			alignment_sh_guru(fastq_map, self._project, self._seq.seq.item(), ref_path, sam_output, ds_sam_output, sh_output)
			logging.info("Alignment jobs are submitted to GURU..")

			# wait for alignment to finish and call mutations

		elif self._env == "BC":
			# make sh files to submit to BC
			sh_output = os.path.join(self._output, "BC_sh")
			os.system("mkdir "+sh_output)

		


def _mut_count(self):
	"""
	Count mutations in input sam files
	"""
	# build lookup table
	pass



def mut_count_sh_guru():
	"""
		"""
	pass


def alignment_sh_guru(fastq_map, ref_name, ref_seq, ref_path, sam_path, ds_sam_path, sh_output):
	"""
	fastq_map: df contains paths to fastq files and downsamled fastq files
	ref_name: name for the reference sequence (same as project name)
	ref_seq: reference sequence to put in fasta file
	ref_path: build fasta in this path
	sam_path: path to sam files
	ds_sam_path: path to downsampled sam files
	sh_output: make sh files to submit to the cluster
	make referencez
	make SGE sh files
	submit sh files to SGE
	monitor jobs
	"""
	## build reference 
	ref = alignment.make_ref(ref_name, ref_seq, ref_path, settings.guru_BOWTIE2_BUILD)

	for index, row in fastq_map.iterrows():
		sample_name = os.path.basename(row["R1"]).split("_")[0]

			shfile = os.path.join(sh_output, sample_name+"_aln.sh") # for each sample, the alignment is for both R1 and R2 (they are aligning separately)
			log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, settings.guru_BOWTIE2, shfile)

			sub_cmd = f"qsub -cwd -N {'aln_'+sample_name} -e {log_file} {shfile}"
			os.system(sub_cmd)

			shfile_ds = os.path.join(sh_output, sample_name+"_aln_ds.sh")
			log_file_ds = alignment.align_main(ref, row["r1_ds"], row["r2_ds"], sam_path, settings.guru_BOWTIE2, shfile_ds)

			sub_cmd = f"qsub -cwd -N {'aln_ds_'+sample_name} -e {log_file_ds} {shfile_ds}"
			os.system(sub_cmd)
			break


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
			break
	# gzip everything 
	os.system(f"gzip {ds_output}/*.fastq")
	fastq_map[["r1_ds", "r2_ds"]] = pd.DataFrame(list(zip(ds_r1files, ds_r2files)))
	return fastq_map


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='TileSeq mutation counts')
	parser.add_argument("-f", "--fastq", help="Path to all fastq files you want to analyze", required= True)
	parser.add_argument("-o", "--output", help="Output folder", required = True)
	parser.add_argument("-log", "--log_level", help="set log level: debug, info, warning, error, critical.", default = "debug")
	parser.add_argument("-p", "--param", help="json paramter file", required = True)
	parser.add_argument("-env", "--environment", help= "The cluster used to run this script", default="GURU")
	parser.add_argument("-n", "--n_reads", help="Used for downsampling the files. n_reads will remain", default = 30000)
	args = parser.parse_args()

	f = args.fastq
	out = args.output
	param = args.param
	env = args.environment
	n = args.n_reads

	log_level = args.log_level

	print(f"Log level: {log_level}")

	# read json file 
	project, seq, cds_seq, tile_map, region_map, samples = help_functions.parse_json(param)

	# Initialize MutCount main 
	mc = MutCount(project, seq, cds_seq, tile_map, region_map, samples, f, out, log_level, env, n)
	# alignment
	# return job ID list (for checking if the jobs are running still)
	job_list = mc._align_sh_()

