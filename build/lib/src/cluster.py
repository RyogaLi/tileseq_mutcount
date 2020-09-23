#!/usr/bin/env python3.6

# The functions in this script is used to manage submitting jobs to different clusters
# This script is called by main.py
import pandas as pd
import os
import re
import subprocess
import time

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
    submit jobs to BC/BC2
    return a df with columns: [R1, R2, r1_sam, r2_sam]
    return a list of job id that we just submited
    """

    # build reference
    ref = alignment.make_ref(ref_name, ref_seq, ref_path, settings.bc2_BOWTIE2_BUILD)
    # store sam paths
    fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
    time = at
    all_job_id = []
    for index, row in fastq_map.iterrows(): # go through all the fastq pairs
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        shfile = os.path.join(sh_output, f"Aln_{sample_name}.sh")
        r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, settings.bc2_BOWTIE2, shfile)
        row["r1_sam"] = r1_sam
        row["r2_sam"] = r2_sam
        # create log file for alignment
        sam_log_f = os.path.join(sam_path, f"{sample_name}.log")
        sub_cmd = ["submitjob2", "-w", str(time), "-c", "4", str(shfile), "2>", sam_log_f]
        jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)

        job_id = jobs.stdout.decode("utf-8").strip().split(".")[0]
        all_job_id.append(job_id)
        logging.info(f"{sample_name}: job id - {job_id}")

    return fastq_map, all_job_id

def alignment_sh_dc(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging):
    """
    submit jobs to BC/BC2
    return a df with columns: [R1, R2, r1_sam, r2_sam]
    return a list of job id that we just submited
    """

    # build reference
    ref = alignment.make_ref(ref_name, ref_seq, ref_path, settings.dc_BOWTIE2_BUILD)
    # store sam paths
    fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
    time = at
    all_job_id = []
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
        job_id = jobs.stdout.decode("utf-8").strip().split(".")[0]
        all_job_id.append(job_id)
        logging.info(f"Sample {sample_name}: job id - {job_id}")

    return fastq_map, all_job_id


def mut_count_sh_bc(sample_name, cmd, mt, sh_output_dir,logger):
    """
    Submit mutation count jobs to BC

    """
    # go through files df and submit jobs for each pair of sam files
    # counting mutations in raw sam output files
    shfile = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.sh")
    with open(shfile, "w") as sh:
        sh.write(cmd+"\n")
        os.system(f"chmod 755 {shfile}")
    # submit this to the cluster
    sub_cmd = ["submitjob2","-w", str(mt), "-c", "1", "-m", "15", shfile]
    job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
    job_id = job.stdout.decode("utf-8").strip().split(".")[0]
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id

def mut_count_sh_dc(sample_name, cmd, mt, sh_output_dir, logger):
    """
    Submit mutation count jobs to DC
    """
    # go through files df and submit jobs for each pair of sam files
    # counting mutations in raw sam output files
    shfile = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.sh")
    with open(shfile, "w") as sh:
        sh.write(cmd+"\n")
        os.system(f"chmod 755 {shfile}")
    #sample_error_file = os.path.join(log_dir, f"sample_{sample_name}.log")
    # submit this to the cluster
    sub_cmd = ["submitjob", "-w", str(mt), "-c", "1", "-m", "15", shfile]
    job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
    job_id = job.stdout.decode("utf-8").strip()
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id


def parse_jobs(job_list, logger):
    """
    return true if all the jobs in job list finished
    else wait for 10 mins and return how man jobs are running and queued
    job_list: list of job ids
    logger: logging object
    """
    qstat_cmd = ["qstat"] + job_list
    job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    qstat_out = job.stdout.decode("utf-8")
    qstat_err = job.stderr.decode("utf-8")

    f_id = []
    updated_list = []
    while True:
        running = []
        queued = []
        completed = []

        if qstat_err != "":
            # jobs might be finished and no longer in queue
            # extract job id from std err
            err = qstat_err.split("\n")[:-1]
            id_regex = re.compile(r"(\d+).bc")
            f_id = [] # finished jobs
            for i in err:
                try:
                    match = id_regex.search(i)
                    job_id = match.group(1)
                    f_id.append(job_id)
                except:
                    print(i)
                    continue
            err_id = set(f_id)
            updated_list = [x for x in job_list if x not in err_id]

        if qstat_out != "":
            qstat_out = qstat_out.split("\n")[:-1]
            id_regex = re.compile(r"(\d+).bc.+(R|Q|C|E)")

            for line in qstat_out:
                if ("---" in line) or ("Job ID" in line): continue
                match = id_regex.search(line)
                job_id = match.group(1)
                job_s = match.group(2)
                if job_s == "E" or job_s == "C":
                    completed.append(job_id)
                elif job_s == "R":
                    running.append(job_id)
                elif job_s == "Q":
                    queued.append(job_id)

        logger.info(f"{len(queued)} jobs queued")
        logger.info(f"{len(running)} jobs running")
        final_list = list(set(updated_list+running+queued))
        if final_list == []:
            return True
        else:
            # check in 10min
            time.sleep(600)
            job_list = final_list
            qstat_cmd = ["qstat"] + job_list
            job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            qstat_out = job.stdout.decode("utf-8")
            qstat_err = job.stderr.decode("utf-8")

if __name__ == "__main__":
    # test job list
    job_list = ["291879", "29171333", "29171340", "29171466"]
    parse_jobs(job_list)
