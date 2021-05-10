#!/usr/bin/env python3.7

# The functions in this script is used to manage submitting jobs to different clusters
# This script is called by main.py
import pandas as pd
import os
import re
import subprocess
import time

# other modules
from TileSeqMut import alignment

# phix reference sequence
phix = ""


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
    ref = alignment.make_ref(ref_name, ref_seq, ref_path)

    for index, row in fastq_map.iterrows():
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        shfile = os.path.join(sh_output, sample_name+"_aln.sh") # for each sample, the alignment is for both R1 and R2 (they are aligning separately)
        log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, shfile)

        sub_cmd = f"qsub -cwd -N {'aln_'+sample_name} -e {log_file} {shfile}"
        os.system(sub_cmd)
        shfile_ds = os.path.join(sh_output, sample_name+"_aln_ds.sh")
        log_file_ds = alignment.align_main(ref, row["r1_ds"], row["r2_ds"], ds_sam_path, shfile_ds)
        sub_cmd = f"qsub -cwd -N {'aln_ds_'+sample_name} -e {log_file_ds} {shfile_ds}"
        os.system(sub_cmd)
        #break

def alignment_sh_galen(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging, rc):
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
    ref = alignment.make_ref(ref_name, ref_seq, ref_path)
    phix = alignment.make_ref("phix", ref_seq, ref_path)

    # store sam paths
    fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
    all_job_id = []
    for index, row in fastq_map.iterrows():
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        # for each sample, the alignment is for both R1 and R2 (they are aligned separately)
        shfile = os.path.join(sh_output, sample_name+"_aln.sh")

        # write header to sh file
        # assume at is single digit
        time_request = f"0{at}:00:00"
        # create log file for alignment
        sam_log_f = os.path.join(sam_path, f"{sample_name}")
        header = f"#!/bin/bash\n#SBATCH --time={time_request}\n#SBATCH --job-name={sample_name}\n#SBATCH " \
                 f"--error={sam_log_f}-%j.log\n#SBATCH --output={sam_log_f}-%j.log\n"

        if "Undetermined" in sample_name:
            r1_sam, r2_sam, log_file = alignment.align_main(phix, row["R1"], row["R2"], sam_path, shfile, rc=rc,
                                                            header=header)
        else:
            r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, shfile, rc=rc, header=header)

        row["r1_sam"] = r1_sam
        row["r2_sam"] = r2_sam
        sub_cmd = ["sbatch", str(shfile)]
        jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = jobs.stdout.decode("utf-8").strip().split(".")[0]
        all_job_id.append(job_id.split()[-1])
        logging.info(f"{sample_name}: job id - {job_id}")
    return fastq_map, all_job_id


def alignment_sh_bc2(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging, rc):
    """
    submit jobs to BC/BC2
    return a df with columns: [R1, R2, r1_sam, r2_sam]
    return a list of job id that we just submited
    """

    # build reference
    ref = alignment.make_ref(ref_name, ref_seq, ref_path)
    # store sam paths
    fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
    time = at
    all_job_id = []
    for index, row in fastq_map.iterrows(): # go through all the fastq pairs
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        shfile = os.path.join(sh_output, f"Aln_{sample_name}.sh")
        r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, shfile, rc=rc)
        row["r1_sam"] = r1_sam
        row["r2_sam"] = r2_sam
        # create log file for alignment
        sam_log_f = os.path.join(sam_path, f"{sample_name}.log")
        sub_cmd = ["submitjob2", "-w", str(time), "-c", "1", str(shfile), "2>", sam_log_f]
        jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)

        job_id = jobs.stdout.decode("utf-8").strip().split(".")[0]
        all_job_id.append(job_id)
        logging.info(f"{sample_name}: job id - {job_id}")

    return fastq_map, all_job_id


def alignment_sh_dc(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging, rc):
    """
    submit jobs to BC/BC2
    return a df with columns: [R1, R2, r1_sam, r2_sam]
    return a list of job id that we just submited
    """

    # build reference
    ref = alignment.make_ref(ref_name, ref_seq, ref_path)
    # store sam paths
    fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
    time = at
    all_job_id = []
    for index, row in fastq_map.iterrows(): # go through all the fastq pairs
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        shfile = os.path.join(sh_output, f"Aln_{sample_name}.sh")
        r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, shfile, rc=rc)
        row["r1_sam"] = r1_sam
        row["r2_sam"] = r2_sam
        # create log file for alignment
        sam_log_f = os.path.join(sam_path, f"{sample_name}.log")
        sub_cmd = ["submitjob","-w", str(time),"-c", "1", str(shfile), "2>", sam_log_f]
        jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = jobs.stdout.decode("utf-8").strip().split(".")[0]
        all_job_id.append(job_id)
        logging.info(f"Sample {sample_name}: job id - {job_id}")

    return fastq_map, all_job_id


def mut_count_sh_bc(sample_name, cmd, mt, mm, sh_output_dir, logger, cores):
    """
    Submit mutation count jobs to BC

    """
    # go through files df and submit jobs for each pair of sam files
    # counting mutations in raw sam output files
    shfile = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.sh")
    log_f = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.log")
    with open(shfile, "w") as sh:
        sh.write(cmd+"\n")
        os.system(f"chmod 755 {shfile}")
    # submit this to the cluster
    sub_cmd = ["submitjob2","-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
    logger.debug(sub_cmd)
    job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
    job_id = job.stdout.decode("utf-8").strip().split(".")[0]
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id


def mut_count_sh_dc(sample_name, cmd, mt, mm, sh_output_dir, logger, cores):
    """
    Submit mutation count jobs to DC
    """
    # go through files df and submit jobs for each pair of sam files
    # counting mutations in raw sam output files
    shfile = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.sh")
    log_f = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.log")
    with open(shfile, "w") as sh:
        sh.write(cmd+"\n")
        os.system(f"chmod 755 {shfile}")
    #sample_error_file = os.path.join(log_dir, f"sample_{sample_name}.log")
    # submit this to the cluster
    sub_cmd = ["submitjob", "-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
    logger.debug(sub_cmd)
    job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
    job_id = job.stdout.decode("utf-8").strip()
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id


def mut_count_sh_galen(sample_name, cmd, mt, mm, sh_output_dir, logger, cores):
    """
    Submit mutation count jobs to DC
    """
    # go through files df and submit jobs for each pair of sam files
    # counting mutations in raw sam output files
    shfile = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.sh")
    log_f = os.path.join(sh_output_dir, f"Mut_count_{sample_name}")
    time_request = f"{mt}:00:00"

    header = f"#!/bin/bash\n#SBATCH --time={time_request}\n#SBATCH --job-name={sample_name}\n#SBATCH " \
             f"--cpus-per-task={cores}\n#SBATCH --error={log_f}-%j.log\n#SBATCH --mem={mm}G\n#SBATCH " \
             f"--output={log_f}-%j.log\n"

    with open(shfile, "w") as sh:
        sh.write(header)
        sh.write(cmd+"\n")
        os.system(f"chmod 755 {shfile}")
    #sample_error_file = os.path.join(log_dir, f"sample_{sample_name}.log")
    # submit this to the cluster
    sub_cmd = ["sbatch", str(shfile)]
    logger.debug(sub_cmd)
    job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
    job_id = job.stdout.decode("utf-8").strip()
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id.split()[-1]


def parse_jobs(job_list, env, logger):
    """
    return true if all the jobs in job list finished
    else wait for 10 mins and return how man jobs are running and queued
    job_list: list of job ids
    logger: logging object
    """
    qstat_cmd = ["qstat"] + job_list
    job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    qstat_out = job.stdout.decode("utf-8", errors="replace")
    qstat_err = job.stderr.decode("utf-8", errors="replace")

    f_id = []
    updated_list = []
    while True:
        running = []
        queued = []
        completed = []

        if qstat_err != "":
            # jobs might be finished and no longer in queue
            # extract job id from std err
            if env == "BC" or env == "BC2":
                err = qstat_err.split("\n")[:-1]
                id_regex = re.compile(r"(\d+).bc")
            elif env == "DC":
                err = qstat_err.split("\n")[:-1]
                id_regex = re.compile(r"(\d+).dc[0-9]+")

            f_id = [] # finished jobs
            for i in err:
                try:
                    match = id_regex.search(i)
                    job_id = match.group(1)
                    f_id.append(job_id)
                except:
                    logger.warning(i)
                    continue
            err_id = set(f_id)
            updated_list = [x for x in job_list if x not in err_id]

        if qstat_out != "":
            qstat_out = qstat_out.split("\n")[:-1]
            if env == "BC" or env == "BC2":
                id_regex = re.compile(r"(\d+).bc.+(R|Q|C|E)")
            elif env == "DC":
                id_regex = re.compile(r"(\d+).dc[0-9]+.+(R|Q|C|E)")

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
            qstat_out = job.stdout.decode("utf-8", errors="replace")
            qstat_err = job.stderr.decode("utf-8", errors="replace")


def parse_jobs_galen(job_list, logger):
    """
    Galen uses slurm scheduler, different from BC and DC
    return true if all the jobs in job list finished
    else wait for 10 mins and return how man jobs are running and queued
    job_list: list of job ids
    logger: logging object
    """
    qstat_cmd = ["squeue", "-j", ",".join(job_list)]
    job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    qstat_out = job.stdout.decode("utf-8", errors="replace")
    qstat_err = job.stderr.decode("utf-8", errors="replace")
    logger.debug(qstat_out)
    logger.debug(qstat_err)
    while True:
        running_jobs, queued_jobs, completed_jobs = [], [], []
        if qstat_out != "":
            # make df
            qstat_df = pd.DataFrame([i.split() for i in qstat_out.split("\n")])
            qstat_df = qstat_df.rename(columns=qstat_df.iloc[0])
            qstat_df = qstat_df.drop(qstat_df.index[1])
            logger.debug(qstat_df)
            # get all active job ID
            running_jobs = qstat_df[qstat_df["ST"] == "R"]["JOBID"].tolist()
            # in queue
            queued_jobs = qstat_df[qstat_df["ST"] == "PD"]["JOBID"].tolist()
            # completing
            completed_jobs = qstat_df[qstat_df["ST"] == "CG"]["JOBID"].tolist()
            logger.debug(running_jobs)
            logger.debug(queued_jobs)

        logger.info(f"{len(queued_jobs)} jobs queued")
        logger.info(f"{len(running_jobs)} jobs running")

        final_list = list(set(running_jobs + queued_jobs))
        if final_list == []:
            return True
        else:
            # check in 10min
            time.sleep(600)
            job_list = final_list
            qstat_cmd = ["squeue", "-j", ",".join(job_list)]
            job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            qstat_out = job.stdout.decode("utf-8", errors="replace")
            qstat_err = job.stderr.decode("utf-8", errors="replace")


def submit_given_jobs(shfile, logger, mt, mm, cores, env=""):
    """
    for a given sh file, submit to cluster, return job id
    """
    sample_name = shfile.split(".")[0].split("_")[-1]
    sh_dir = os.path.dirname(shfile)
    log_f = os.path.join(sh_dir, f"Mut_count_{sample_name}.log")
    if env == "GALEN":
        sub_cmd = ["sbatch", str(shfile)]
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip().split()[-1]
    elif env == "BC2" or env == "BC":
        sub_cmd = ["submitjob2", "-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip()
    elif env == "DC":
        sub_cmd = ["submitjob", "-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip()
    else:
        raise ValueError("Wrong environment code")
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id


if __name__ == "__main__":
    # test job list
    job_list = ["42091", "42090", "42089", "42080"]
    # parse_jobs(job_list, "DC", "")
    parse_jobs_galen(job_list, "")
