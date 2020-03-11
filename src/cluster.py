#~/lib/Python-3.6.4/python

# The functions in this script is used to manage submitting jobs to different clusters
# This script is called by main.py
import pandas as pd
import os
import subprocess

# other modules
import alignment
import help_functions
import settings

def alignment_sh_guru(fastq_map, ref_name, ref_seq, ref_path, sam_path, ds_sam_path, sh_output):
	"""
	fastq_map: df contains paths to fastq files and downsamled fastq files
	ref_name: name for the reference sequence (same as project name)
	ref_seq: reference sequence to put in fasta file
	ref_path: build fasta in this path
	sam_path: path to sam files
	ds_sam_path: path to downsampled sam files
	sh_output: make sh files to submit to the cluster

	make reference for bowtie2
	make SGE sh files
	submit sh files to SGE
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
		log_file_ds = alignment.align_main(ref, row["r1_ds"], row["r2_ds"], ds_sam_path, settings.guru_BOWTIE2, shfile_ds)
		sub_cmd = f"qsub -cwd -N {'aln_ds_'+sample_name} -e {log_file_ds} {shfile_ds}"
		os.system(sub_cmd)
		#break

def alignment_sh_bc2(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging):
	"""
	"""

	# build reference
	ref = alignment.make_ref(ref_name, ref_seq, ref_path, settings.bc2_BOWTIE2_BUILD)
	# store sam paths
	fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
	time = at
	for index, row in fastq_map.iterrows(): # go through all the fastq pairs
		sample_name = os.path.basename(row["R1"]).split("_")[0]

		shfile = os.path.join(sh_output, f"Aln_{sample_name}.sh")
		r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, settings.bc2_BOWTIE2, shfile)
		row["r1_sam"] = r1_sam
		row["r2_sam"] = r2_sam
		# create log file for alignment
		sam_log_f = os.path.join(sam_path, f"{sample_name}.log")
		sub_cmd = ["submitjob", str(time), "-c", "4", str(shfile), "2>", sam_log_f]
		jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
		ids = jobs.stdout.decode("utf-8").strip()
		logging.info(f"{sample_name}: job id - {ids}")

	return fastq_map

def alignment_sh_dc(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging):
	"""
		"""

	# build reference
	ref = alignment.make_ref(ref_name, ref_seq, ref_path, settings.dc_BOWTIE2_BUILD)
	# store sam paths
	fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
	time = at
	for index, row in fastq_map.iterrows(): # go through all the fastq pairs
		sample_name = os.path.basename(row["R1"]).split("_")[0]

		shfile = os.path.join(sh_output, f"Aln_{sample_name}.sh")
		r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, settings.dc_BOWTIE2, shfile)
		row["r1_sam"] = r1_sam
		row["r2_sam"] = r2_sam
		# create log file for alignment
		sam_log_f = os.path.join(sam_path, f"{sample_name}.log")
		sub_cmd = ["submitjob","-w", str(time),"-c", "4", str(shfile), "2>", sam_log_f]
		jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
		ids = jobs.stdout.decode("utf-8").strip()
		logging.info(f"Sample {sample_name}: job id - {ids}")

	return fastq_map


def mut_count_sh_bc(files_df, output_dir, param_json, sh_output, min_cover, mt, log_dir, logging, qual):
	"""
	Submit mutation count jobs to BC
	files_df: dataframe contains full path to all the sam files
	output_dir: path to save the mutation count output
	param_json: json parameter file
	sh_output: path to save sh files (you can find all the executable bash scripts for all the samples)
	min_cover: minimum % coverage of the tile
	mt: time needed for this run
	log_dir: directory to save the mutation count log file
	logging: main logging object
	qual: quality filter (posterior prob cut off)
	"""
	# go through files df and submit jobs for each pair of sam files
	py_path = os.path.abspath("count_mut.py")
	time = mt
	for index, row in files_df.iterrows():
			sample_name = os.path.basename(row["r1_sam"]).split("_")[0]
			# counting mutations in raw sam output files
			shfile = os.path.join(sh_output, f"Mut_count_{sample_name}.sh")
			cmd = f"python {py_path} -r1 {row['r1_sam']} -r2 {row['r2_sam']} -o {output_dir} -p {param_json} -qual {qual} -mutlog {log_dir} -min {min_cover}"
			with open(shfile, "w") as sh:
					sh.write(cmd+"\n")
					os.system(f"chmod 755 {shfile}")
			#sample_error_file = os.path.join(log_dir, f"sample_{sample_name}.log")
			# submit this to the cluster
			sub_cmd = ["submitjob", str(time), "-c", "8", "-m", "4", shfile]
			job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
			ids = job.stdout.decode("utf-8").strip()
			# log sample name and job id
			logging.info(f"Sample {sample_name}: job id - {ids}")

def mut_count_sh_dc(files_df, output_dir, param_json, sh_output, min_cover, mt, log_dir, logging, qual):
	"""
	Submit mutation count jobs to BC
	files_df: dataframe contains full path to all the sam files
	output_dir: path to save the mutation count output
	param_json: json parameter file
	sh_output: path to save sh files (you can find all the executable bash scripts for all the samples)
	min_cover: minimum % coverage of the tile
	mt: time needed for this run
	log_dir: directory to save the mutation count log file
	logging: main logging object
	qual: quality filter (posterior prob cut off)
	"""
	# go through files df and submit jobs for each pair of sam files
	py_path = os.path.abspath("count_mut.py")
	time = mt
	for index, row in files_df.iterrows():
			sample_name = os.path.basename(row["r1_sam"]).split("_")[0]
			# counting mutations in raw sam output files
			shfile = os.path.join(sh_output, f"Mut_count_{sample_name}.sh")
			cmd = f"python {py_path} -r1 {row['r1_sam']} -r2 {row['r2_sam']} -o {output_dir} -p {param_json} -qual {qual} -mutlog {log_dir} -min {min_cover}"
			with open(shfile, "w") as sh:
					sh.write(cmd+"\n")
					os.system(f"chmod 755 {shfile}")
					#sample_error_file = os.path.join(log_dir, f"sample_{sample_name}.log")
					# submit this to the cluster
					sub_cmd = ["submitjob", "-w", str(time),"-c", "8", "-m", "4", shfile]
					job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
					ids = job.stdout.decode("utf-8").strip()
					# log sample name and job id
					logging.info(f"Sample {sample_name}: job id - {ids}")

