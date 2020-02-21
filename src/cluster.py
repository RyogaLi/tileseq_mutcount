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

def alignment_sh_bc2(fastq_map, ref_name, ref_seq, ref_path, sam_path, ds_sam_path, sh_output, logging):
	"""
	"""
	# build reference
	ref = alignment.make_ref(ref_name, ref_seq, ref_path, settings.bc2_BOWTIE2_BUILD)
	# store sam paths
	fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam", "r1_sam_ds", "r2_sam_ds"])])
	for index, row in fastq_map.iterrows(): # go through all the fastq pairs
		sample_name = os.path.basename(row["R1"]).split("_")[0]

		shfile = os.path.join(sh_output, f"Aln_{sample_name}.sh")
		r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, settings.bc2_BOWTIE2, shfile)
		row["r1_sam"] = r1_sam
		row["r2_sam"] = r2_sam
		# create log file for alignment
		sam_log_f = os.path.join(sam_path, f"{sample_name}.log")
		time = 8 # schedule this alignment for 8 hours (this is more than what we need)
		sub_cmd = ["submitjob",str(time), str(shfile), "2>", sam_log_f]
		jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
		ids = jobs.stdout.decode("utf-8").strip()
		logging.info(f"{ids}")
		time_ds = 1
		shfile_ds = os.path.join(sh_output, f"Aln_ds_{sample_name}.sh")
		r1_sam_ds, r2_sam_ds, log_file_ds = alignment.align_main(ref, row["r1_ds"], row["r2_ds"], ds_sam_path, settings.bc2_BOWTIE2, shfile_ds)
		row["r1_sam_ds"] = r1_sam_ds
		row["r2_sam_ds"] = r2_sam_ds
		# create log file for alignment (ds)
		sam_log_ds = os.path.join(ds_sam_path, f"{sample_name}_ds.log")
		sub_cmd_ds = ["submitjob", str(time_ds), str(shfile_ds), "2>", sam_log_ds]
		jobs_ds = subprocess.run(sub_cmd_ds, stdout=subprocess.PIPE)
		id_ds = jobs_ds.stdout.decode("utf-8").strip()
		logging.info(f"{id_ds}")
		#break
	# return the merged df
	return fastq_map


def mut_count_sh_bc(files_df, output_dir, param_json, sh_output, log_dir, logging):
	"""
	Submit mutation count jobs to BC
	files_df: dataframe contains full path to all the sam files
	output_dir: path to save the mutation count output
	param_json: json parameter file
	sh_output: path to save sh files (you can find all the executable bash scripts for all the samples)
	log_dir: directory to save the mutation count log file
	logging: main logging object
	"""
	# go through files df and submit jobs for each pair of sam files
	py_path = os.path.abspath("count_mut.py")
	for index, row in files_df.iterrows():
		sample_name = os.path.basename(row["r1_sam"]).split("_")[0]
		# create log files
		log_f = os.path.join(log_dir, "sample_"+str(sample_name)+"_mut.log")
		#log_f_ds = os.path.join(log_dir, "sample_"+str(sample_name)+"_ds_mut.log")

		# counting mutations in raw sam output files
		time = 20
		shfile = os.path.join(sh_output, f"Mut_count_{sample_name}.sh")
		cmd = f"python {py_path} -r1 {row['r1_sam']} -r2 {row['r2_sam']} -o {output_dir} -p {param_json} -logf {log_f}"
		with open(shfile, "w") as sh:
			sh.write(cmd+"\n")
		os.system(f"chmod 755 {shfile}")

		# submit this to the cluster
		sub_cmd = ["submitjob", str(time), "-m", "10", shfile, "2>>", log_f]
		job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
		ids = job.stdout.decode("utf-8").strip()
		# log sample name and job id
		logging.info(f"{sample_name}: {ids}")
		print(f"{sample_name}: {ids}")
		# counting mutations in downsampled sam output files
		#time = 1
		#shfile_ds = os.path.join(sh_output, f"Mut_count_{sample_name}_ds.sh")
		#cmd = f"python {py_path} -r1 {row['r1_sam_ds']} -r2 {row['r2_sam_ds']} -o {output_dir} -p {param_json} -logf {log_f_ds}"
		#with open(shfile_ds, "w") as sh_ds:
		#	sh_ds.write(cmd+"\n")
		#os.system(f"chmod 750 {shfile_ds}")

		#sub_cmd_ds = ["submitjob", str(time), shfile_ds, "2>>", log_f]
		#jobs_ds = subprocess.run(sub_cmd_ds, stdout=subprocess.PIPE)
		#id_ds = jobs_ds.stdout.decode("utf-8").strip()
		# log sample name and job ID
		#logging.info(f"{id_ds}")

