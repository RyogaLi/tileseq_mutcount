#!/usr/bin/env python3.7

# Main script for sequencing analysis
# Input of this script:
#        - JSON file containing input parameters
#    - Fastq files (gzipped)

# what does this script do?
# 1. Read input csv file and convert to json (using csv2Json.R (you need to
# install required packages in R)
# 1. Read paramters from json file
# 2. Read input fastq file, align reads to reference sequence, generate sam files. (bowtie2 required)
# 4. From sam files, count mutations
# 5. Output mutation counts to summary.csv

# modules
import re
import subprocess

import pandas as pd
import numpy as np
import os
import glob
import argparse
import logging
import datetime
import shutil

# pakage modules
from TileSeqMut import help_functions
from TileSeqMut import count_mut
from TileSeqMut import cluster


class fastq2counts(object):
    """Main class to read fastq files and output mutation counts."""

    def __init__(self, param_path, output_dir, main_log, args):
        """
        Initialize mutation counts
        param_path: Full path to param.json file
        output_dir: Directory to save all the output files
        main_log: logging object
        args: user inpu arguments
        """
        # path to this python file
        self._main_path = os.path.abspath(__file__)
        # logger
        self._logging = main_log
        # parameter json file
        self._param_json = param_path
        # user input arguments
        self._args = args
        # convert environment input into upper case
        self._args.environment = self._args.environment.upper()

        # if input arguments provided fastq path, get all fastq files
        if self._args.fastq:
            self._fastq_list = glob.glob(self._args.fastq+"/*.fastq.gz")
        else:
            self._fastq_list = []
        # output dir for this run
        self._output = output_dir

        # parse parameter json file
        self._project, self._seq, self._cds_seq, self._tile_map, self._region_map, \
            self._samples, self._var, self._relations = help_functions.parse_json(param_path)
        self._sample_names = self._samples["Sample ID"].tolist()
        self._project = self._project.replace(" ", "_")  # project name
        # make main log
        self._log = self._logging.getLogger("main.log")
        # self._log.debug(self._args.environment)

    def _init_skip(self, skip, r1=None, r2=None):
        """
        if skip == T: skip alignment
            make suboutput dir with time stamp
        else: not skip alignment
            make main output dir with time stamp
        """
        # initialize read 1/read 2
        self._r1 = ""
        self._r2 = ""
        self._skip = False
        if skip:  # only run mutation counts on the aligned samples
            self._skip = True
            if self._args.r1 and self._args.r2:
                self._r1 = self._args.r1
                self._r2 = self._args.r2
            else:
                self._log.info(f"Sam files are read from {self._output}")

    def _align_sh_(self):
        """
        1. Align each fastq file to reference (submit one job for each alignment)
        Output sam files into output folders (sam_files)
        """
        # VALIDATE: if all samples are present in fastq files
        # also check if any fastq file is empty
        fastq_sample_id = [os.path.basename(i).split("_")[0] for i in self._fastq_list]
        fastq_sample_id = list(set(fastq_sample_id))
        self._log.debug(f"Fastq sample ID in input folder: {fastq_sample_id}")
        # validation
        # check if input fastq files contains all the samples in paramter file
        self._log.debug(f"Sample names in parameter sheet: {set(self._sample_names)}")
        if set(self._sample_names).issubset(fastq_sample_id):
            self._log.info("All samples found")
            self._log.info(f"In total there are {len(list(set(self._sample_names)))} samples in the csv file")
            self._log.info(f"In total there are {len(fastq_sample_id)} fastq files")
        else:
            test= list(np.setdiff1d(self._sample_names,fastq_sample_id))
            join_list = ",".join(test)
            self._log.error("fastq files do not match input samples.")
            self._log.error(f"Fastq files not found for sample {join_list}")
            raise FileNotFoundError()

        # create mapping for R1 and R2 for each sample
        # ONLY samples provided in the parameter files will be analyzed
        fastq_map = []
        for r1 in self._fastq_list:
            ID = os.path.basename(r1).split("_")[0]
            if ("_R1_" in r1) and (ID in self._sample_names):
                    r2 = r1.replace("_R1_", "_R2_")
                    fastq_map.append([r1,r2])

        # make folder to store alignment sam files
        sam_output = os.path.join(self._output, "sam_files/")
        os.system("mkdir " + sam_output)

        # make folder to store ref
        ref_path = os.path.join(self._output, "ref/")
        os.system("mkdir " + ref_path)

        # convert to df
        fastq_map = pd.DataFrame(fastq_map, columns=["R1", "R2"])
        # if phix
        phix_fasta = os.path.join(os.path.dirname(self._main_path), "data/phix.fasta")
        # self._log.info(os.path.dirname(self._main_path))
        # save this to ref path
        # self._log.info(phix_fasta)
        cmd = f"cp {phix_fasta} {ref_path}"
        # self._log.info(cmd)
        os.system(cmd)
        self._phix_fasta = []
        if self._args.calibratePhredPhix:
            # check if Undetermined reads are in the same folder
            self._phix_fastq = glob.glob(self._args.fastq+"/Undetermined_*.fastq.gz")
            if len(self._phix_fastq) < 2:
                raise ValueError("Cannot find Undetermined fastq files for phix alignment")
        #if self._phix_fastq != []:
            # add this to the last row
            fastq_map.loc[len(fastq_map)] = self._phix_fastq

        rc = False
        if self._args.rc:
            rc = True

        # GURU # depreciated
        # if self._args.environment == "GURU":
        #     # make sh files to submit to GURU
        #     sh_output = os.path.join(self._output, "GURU_jobs")
        #     os.system("mkdir "+sh_output)
        #     # for each pair of fastq files, submit jobs to run alignment
        #     # it takes the following arguments: R1, R2, ref, Output sam
        #     # the output log (bowtie log) would be in the same dir
        #     logging.info("Writing sh files for alignment (GURU)")
        #     cluster.alignment_sh_guru(fastq_map, self._project, self._seq.seq.item(), ref_path, sam_output, sh_output, logging)
        #     logging.info("Alignment jobs are submitted to GURU..")
        #     # wait for alignment to finish and call mutations

        if self._args.environment == "GALEN":
            # make sh files to submit to GURU
            sh_output = os.path.join(self._output, "GALEN_jobs")
            os.system("mkdir "+sh_output)
            # for each pair of fastq files, submit jobs to run alignment
            # it takes the following arguments: R1, R2, ref, Output sam
            # the output log (bowtie log) would be in the same dir
            logging.info("Writing sh files for alignment (GALEN)")
            sam_df, job_list = cluster.alignment_sh_galen(fastq_map, self._project, self._seq.seq.values.item(), ref_path,
                                            sam_output, sh_output, self._args.at, self._log, rc)
            logging.info("Alignment jobs are submitted to GALEN..")

        elif self._args.environment == "BC2" or self._args.environment == "BC" or self._args.environment == "DC":
            # make sh files to submit to BC
            sh_output = os.path.join(self._output, "CCBR_sh")
            os.system("mkdir " + sh_output)
            # make sh files to submit to BC
            self._log.info("Submitting alignment jobs...")
            # alignment_sh_bc2(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging)

            sam_df, job_list = cluster.alignment_sh_ccbr(fastq_map, self._project, self._seq.seq.values.item(),
                                                        ref_path, sam_output, sh_output, self._args.at,
                                                        self._log, rc, self._args.environment)
            self._log.info("Alignment jobs are submitte to BC2. Check pbs-output for STDOUT/STDERR")

        else:
            raise ValueError("Wrong environment name")

        self._log.info(f"Total jobs submitted: {len(job_list)}")
        if self._args.environment == 'GALEN':
            finished = cluster.parse_jobs_galen(job_list, self._logging.getLogger("track.jobs"))
        else:
            finished = cluster.parse_jobs(job_list, self._args.environment, self._logging.getLogger("track.jobs"))  #

        if finished:
            self._log.info(f"Alignment jobs are finished!")
            self._log.info(f"Merging alignment log files...")
            # go through alignment output dir and parse all the log files
            log_files = os.path.join(sam_output, "*.log")
            log_list = glob.glob(log_files)
            alignment_master_log = []
            for f in log_list:
                f_info = []
                sample_id = os.path.basename(f).split(".")[0]
                f_info.append(sample_id)
                with open(f, "r") as fp:
                    for i, line in enumerate(fp):
                        # pick corresponding lines in the alignment log file
                        if i == 0 or i == 5 or i == 6 or i == 11:
                            line = line.split(" ")
                            f_info.append(line[0])
                alignment_master_log.append(f_info)
            df = pd.DataFrame(alignment_master_log)
            df.columns = ["sample_id", "R1_reads", "R1_alignment_rate", "R2_reads", "R2_alignment_rate"]
            alignment_log = os.path.join(sam_output, "aligned_rate.log")
            df.to_csv(alignment_log, index=False)
            self._log.info(f"Alignment log file created")

        return sam_output

    def _mut_count(self):
        """
        Count mutations in input sam files
        1. Log path to sam files and downsampled sam files
        2. If sam df not provided, make df for all the sam files
        3. For each pair of sam files, submit jobs to the cluster
        4. Log to main log
        """
        self._log.info("Counting mutations ...")

        # submit job with main.py -r1 and -r2
        # run main.py with -r1 and -r2
        logger_mut = self._logging.getLogger("count.mut.log")
        self._log.info("Removing existing phred files...")
        #phred_flist = glob.glob(os.path.join(self._output, '*_phred.*'))
        #if phred_flist != []:
        #    os.system(f"rm {os.path.join(self._output, '*_phred.*')}")
        mut_counts = count_mut.readSam(self._r1, self._r2, self._param_json, self._args, self._output, self._args.c, logger_mut)
        self._log.info("Running multi-core analysis ... ")
        if self._args.calibratePhredWT:
            self._log.info("Adjusting phred scores using WT samples...")
            adjusted_phred = mut_counts.adjust_er(wt_override=self._args.wt_override)
        elif self._args.calibratePhredPhix:
            self._log.info("Adjusting phred scores using Phix...")
            adjusted_phred = mut_counts.adjust_er_phix()
        else:
            adjusted_phred = []
        mut_counts.multi_core(adjusted_phred)

    def _makejobs(self, sh_output, sam_dir, samples):
        """
        For each pair of sam files in output/sam_files/
        submit mut count job to the cluster
        """
        if sam_dir == "":
            # get samples in param file
            sam_dir = os.path.join(self._args.output, "sam_files/") # read sam file from sam_file
        # get sam files from parameter file
        if not os.path.isdir(sam_dir):
            self._log.error(f"Directory: ./sam_files/ not found in {self._output}")
            raise ValueError()

        self._log.debug(f"Sam files are read from {sam_dir}")
        job_list = []
        for sample in samples:
            # assume all the sam files have the same name format (id_*.sam)
            sam_f_r1 = glob.glob(f"{sam_dir}/{sample}_*R1_*.sam")
            sam_f_r2 = glob.glob(f"{sam_dir}/{sample}_*R2_*.sam")
            if len(sam_f_r1) == 0 or len(sam_f_r2) == 0:
                self._log.error(f"SAM file for sample {sample} not found. Please check your parameter file")
                raise ValueError()

            else:
                self._log.info(f"Sample {sample}")
                self._log.info(f"Read1: {sam_f_r1[0]}")
                self._log.info(f"Read2: {sam_f_r2[0]}")
                sam_id_1 = os.path.basename(sam_f_r1[0]).split("_")[0]
                sam_id_2 = os.path.basename(sam_f_r2[0]).split("_")[0]
                if (sam_id_1 != sample) or (sam_id_1 != sam_id_2) or (sam_id_2 != sample):
                    self._log.error("IDs in sam files don't match!")
                self._r1 = sam_f_r1[0]
                self._r2 = sam_f_r2[0]
            
            job_id = self._submit_mutcount_jobs(sample, sh_output)
            job_list.append(job_id)
        self._log.info("all samples submitted")

        # after submitting all jobs, track all submitted jobs on cluster
        jobs = ",".join(job_list)
        self._log.debug(f"All jobs: {jobs}")
        self._log.info(f"Total jobs running: {len(job_list)}")
        if self._args.environment == 'GALEN':
            finished = cluster.parse_jobs_galen(job_list, self._logging.getLogger("track.jobs"))
        else:
            finished = cluster.parse_jobs(job_list, self._args.environment, self._logging.getLogger("track.jobs"))

        return finished

    def _submit_mutcount_jobs(self, sample, sh_output):
        """
        Make mutation counting jobs and submit for this sample
        """
        # submit job with main.py -r1 and -r2
        # run main.py with -r1 and -r2
        cmd = f"tileseq_mut -n {self._args.name} -r1 {self._r1} -r2 {self._r2} -o {self._output} -p" \
                  f" {self._param_json} --skip_alignment -log {self._args.log_level} -env {self._args.environment} -at {self._args.at} -mt {self._args.mt} -c {self._args.c} " 
        if self._args.sr_Override:
            cmd = cmd + "-override "

        if self._args.wt_override:
            cmd = cmd + "--wt_override "
        
        if self._args.calibratePhredPhix:
            cmd = cmd + "--calibratePhredPhix "

        if self._args.calibratePhredWT:
            cmd = cmd + "--calibratePhredWT "

        if self._args.environment == "BC2" or self._args.environment == "BC" or self._args.environment == "DC":
            logging.info("Submitting mutation counts jobs....")
            job_id = cluster.mut_count_sh_ccbr(sample, cmd, self._args.mt, self._args.mm, sh_output, self._log,
                                             self._args.c, self._args.environment)
        # elif :
        #     logging.info("Submitting mutation counts jobs to DC...")
        #     # (sample_name, cmd, mt, sh_output_dir, logger)
        #     job_id = cluster.mut_count_sh_dc(sample, cmd, self._args.mt, self._args.mm, sh_output, self._log,
        #                                      self._args.c)
        elif self._args.environment == "GALEN":
            logging.info("Submitting mutation counts jobs to GALEN...")
            # (sample_name, cmd, mt, sh_output_dir, logger)
            job_id = cluster.mut_count_sh_galen(sample, cmd, self._args.mt, self._args.mm, sh_output, self._log,
                                                self._args.c)  # this
        else:
            raise ValueError("Wrong environment")

        return job_id

    def _checkoutput(self, f):
        """
        Check how many variants are generated in each file
        """
        mut_n = 0
        with open(f, "r") as mut_output:
            for line in mut_output:
                # skip header
                if "c." in line:
                    mut_n += 1
        return mut_n

    def main(self):
        """
        Main for the fastq2counts object
        """
        sam_dir = ""
        if self._skip == False: # submit jobs for alignment
            # get sam_df from self._align_sh()
            sam_dir = self._align_sh_()

        if self._r1 != "" and self._r2 != "":
            # call functions in count_mut.py
            # init(sam_r1, sam_r2, param, args, ouptut_dir, logger)
            self._mut_count()

        # submit jobs for mutation counting
        # if user did not provide r1 and r2 SAM file
        # get r1 and r2 sam files from output dir provided by user
        elif self._r1 == "" and self._r2 == "":
            # resubmit is used when some jobs in the queue did not finish properly
            if self._args.resubmit:
                # resubmit failed jobs in existing mut_count dir
                # self.output is the mut_count dir in this case
                # find out which jobs failed by going through all counts files

                # first find if any phred calibration files are empty
                # remove empty files
                self._log.info("Removing existing phred files if they are empty...")
                phred_flist = glob.glob(os.path.join(self._output, '*_phred.*'))
                for ph in phred_flist:
                    if os.stat(ph).st_size == 0:
                        os.system(f"rm {os.path.join(self._output, ph)}")
                # get all the mut count files
                mutcount_list = glob.glob(os.path.join(self._output, "counts_sample_*.csv"))
                failed_samples = []
                all_samples_with_csv = []
                # go through all the mut count files
                for f in mutcount_list:
                    m = re.search('.*counts_sample_(.+?).csv', f)
                    # find sample ID in file name
                    if m:
                        sample_id = m.group(1)
                    else:
                        self._log.error("Sample ID not found")
                        raise ValueError(f"Sample ID not found for {f}")
                    all_samples_with_csv.append(sample_id)
                    # check if the file is empty
                    mut_n = self._checkoutput(f)
                    if mut_n == 0:
                        failed_samples.append(sample_id)
                # if any sample does not have a csv file generated (jobs were killed before started)
                other_missing = [i for i in self._sample_names if i not in all_samples_with_csv]
                failed_samples += other_missing
                self._log.info(f"Failed samples: {failed_samples}")
                self._log.info("Resubmitting samples ...")
                # find the path to jobs dir
                env_jobs_dir = os.path.join(self._output, f"{self._args.environment}_jobs")
                parent_dir = os.path.abspath(os.path.join(self._output, os.pardir))
                sam_dir = os.path.join(parent_dir, "sam_files")
                # remake sh files for these samples and submit them
                finished = self._makejobs(env_jobs_dir, sam_dir, failed_samples)
                
                self._log.info("ALL DONE!")
                return finished

            # start alignment jobs
            if self._skip == False:
                time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
                # make mut count dir alignemnt not skipped
                # otherwise the user input dir should be the output dir generated
                self._output = os.path.join(self._output, self._args.name + "_" + time_now + "_mut_count")
                os.makedirs(self._output)
                # align to phix if specified
                if self._args.calibratePhredPhix:
                    # check if Undetermined reads are in the same folder
                    self._phix_fastq = glob.glob(self._args.fastq+"/Undetermined_*.fastq.gz")
                    if len(self._phix_fastq) < 2:
                        raise ValueError("Cannot find Undetermined fastq files for phix alignment")

            # output directory is the mut_count dir
            # make folder to store all the sh files as well as all the log files
            if self._args.environment == "BC2" or self._args.environment == "BC":
                sh_output = os.path.join(self._output, "BC_jobs")

            elif self._args.environment == "DC":
                sh_output = os.path.join(self._output, "DC_jobs")

            elif self._args.environment == "GALEN":
                sh_output = os.path.join(self._output, "GALEN_jobs")

            else:
                self._log.error("Please provide valid environment: BC/BC2/DC/Galen")
                raise ValueError("Please provide valid environment: BC/BC2/DC/Galen")

            self._log.info(f"Mutation count sh files are made in {sh_output}")

            os.mkdir(sh_output)
            finished = self._makejobs(sh_output, sam_dir, self._sample_names)

            if finished:
                self._log.info(f"Mutation counting jobs are finished!")
                self._log.info("Check mutation counts file ...")
                mutcount_list = glob.glob(os.path.join(self._output, "counts_sample_*.csv"))
                self._log.info(f"{len(mutcount_list)} mutation counts file generated")
                if len(self._sample_names) > len(mutcount_list):
                    self._log.error("Job finished but some variant call files are missing!!")
                    missing = list(set(self._sample_names)-set(mutcount_list))
                    missing = ",".join(missing)
                    self._log.error(f"Missing samples: {missing}")
                # delay here for the system to write to the files
                os.system("sleep 180")
                for f in mutcount_list:
                    mut_n = self._checkoutput(f)
                    if mut_n == 0:
                        self._log.error(f"{f} has 0 variants! Check mut log for this sample.")
                    else:
                        self._log.info(f"{f} has {mut_n} variants")


def check(args):
    """
    Validate args and convert csv to JSON
    return path to json
    """

    # two script needed: csv2json.R and calibratePhred.R
    test_csv2json = subprocess.getstatusoutput("csv2json.R -h")
    test_caliPhred = subprocess.getstatusoutput("calibratePhred.R -h")
    if test_csv2json[0] != 0:
        raise OSError(f"csv2json.R failed with error {test_csv2json[1]}")
    if test_caliPhred[0] != 0:
        raise OSError(f"calibratePhred.R failed with error {test_caliPhred[1]}")

    # check bowtie2 and bowtie2-build
    test_bt = subprocess.getstatusoutput("bowtie2 -h")
    test_bt_build = subprocess.getstatusoutput("bowtie2-build -h")
    if test_bt[0] != 0:
        raise OSError(f"bowtie2 failed with error {test_bt[1]}")
    if test_bt_build[0] != 0:
        raise OSError(f"bowtie2-build failed with error {test_bt_build[1]}")

    # check if the correct args are provided
    # if you don't specify --skip_alignment then you cannot provide r1 and r2 sam_files
    if not args.skip_alignment:
        if args.r1 or args.r2:
            raise ValueError(f"Invalid paramters! Please specify --skip_alignment if you want to analyze\
            one pair of sam files")

    # check if output dir exists
    if not os.path.isdir(args.output):
        raise FileNotFoundError(f"Output directory not found: {args.output}")

    # check if fastq dir exists
    if args.fastq and not os.path.isdir(args.fastq):
        raise FileNotFoundError(f"Fastq file path not found: {args.fastq}")

    # try convert csv to json in the same dir as the csv file
    # convert csv file to json
    if args.param.endswith(".csv"):
        print("parameter sheet in csv format, converting to json")
        param_json = args.param.replace(".csv", ".json")
        if os.path.isfile(param_json):
             os.remove(param_json)
        if args.sr_Override:
            convert = f"csv2json.R {args.param} -o {param_json} --srOverride"
        else:
            convert = f"csv2json.R {args.param} -o {param_json}"
        os.system(convert)

    # if the file ends with .json, do nothing
    elif args.param.endswith(".json"):
        param_json = args.param
    else:
        raise ValueError("Please provide valid paramter file format (csv or json)")

    if not os.path.isfile(param_json):
        raise ValueError("Json file does not exist, check conversion!")

    # check environment input
    if args.environment.upper() not in ["BC", "BC2", "DC", "GALEN"]:
        raise ValueError(f"Wrong input for environment: {args.environment}")

    return param_json


def main(args, v):
    """
    Main for fastq2counts
    """
    # get time stamp for this object
    time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    # check input arguments
    param_json = check(args)

    # get basename for the param file
    param_base = os.path.basename(param_json)
    if args.skip_alignment: # if we want to skip the alignment part
        # if user provided R1sam + R2sam
        if args.r1 and args.r2:
            # analyze one pair of r1 and r2 file
            # VALIDATE if both files are sam files
            if not args.r1.endswith(".sam") or not args.r2.endswith(".sam"):
                raise ValueError("Please provide SAM files")
            # copy the parameter file to output dir if the json file is not in the dir
            param_path = os.path.join(args.output, param_base)
            if not os.path.isfile(param_path):
                param_path = shutil.copy(param_json, args.output, follow_symlinks=True)
            # set up loggingthis creates main.log in the mut_count output dir
            main_log_f = os.path.join(args.output, "mutcount_main.log")
            main_log = help_functions.logginginit(args.log_level, main_log_f)

            # initialize mutcount object
            mc = fastq2counts(param_path, args.output, main_log, args)
            mc._init_skip(skip=True, r1=args.r1, r2=args.r2)

        else:
            if not args.resubmit:
                # user provided a csv file and path to sam files
                # for each pair of sam file we would submit one job to the cluster for mutation counting
                # make time stamped output folder for this project
                updated_out = os.path.join(args.output, args.name + "_" + time_now + "_mut_count")
                # sam_dir = os.path.join(args.output, "sam_files/") # read sam file from sam_file
                os.makedirs(updated_out)
            else:
                updated_out = args.output

            args_log_path = os.path.join(updated_out, f"{time_now}_args.log")
            # record input parameters
            write_param(args_log_path, args, v)

            # load json file
            param_path = os.path.join(updated_out, param_base)
            # copy this file to the output folder
            if not os.path.isfile(param_path):
                param_path = shutil.copy(param_json, updated_out, follow_symlinks=True)
            main_log_f = os.path.join(args.output, "mutcount_main.log")
            main_log = help_functions.logginginit(args.log_level, main_log_f)
            # initialize mutcount object
            mc = fastq2counts(param_path, updated_out, main_log, args)
            # skip = skip alignment
            mc._init_skip(skip=True)
    else:
        # alignment
        # make time stamped output folder for this project
        updated_out = os.path.join(args.output, args.name + "_" + time_now)
        os.makedirs(updated_out)  # make directory to save this run
        # write args to file
        args_log_path = os.path.join(updated_out, f"{time_now}_args.log")
        write_param(args_log_path, args, v)

        param_path = os.path.join(updated_out, param_base)
        if not os.path.isfile(param_path):
            param_path = shutil.copy(param_json, updated_out, follow_symlinks=True)

        main_log_f = os.path.join(updated_out, "main.log")
        main_log = help_functions.logginginit(args.log_level, main_log_f)

        # initialize mutcount object
        mc = fastq2counts(param_path, updated_out, main_log, args)
        mc._init_skip(skip=False)

    # start the run
    mc.main()


def write_param(args_log_path, args, v):
    """
    Write input arguments to param.log
    """
    with open(args_log_path, "w") as args_log:
        args_log.write(f"{v}\n")
        for arg in vars(args):
            args_log.write(arg+",")
            args_log.write(f"{getattr(args, arg, 'N/A')}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='TileSeq mutation counts')
    # user input arguments
    parser.add_argument("-f", "--fastq", help="Path to all fastq files you want to analyze", type=str)
    parser.add_argument("-o", "--output", help="Output folder", type=str, required=True)
    parser.add_argument("-p", "--param", help="csv paramter file", type=str, required=True)
    parser.add_argument("-n", "--name", help="Name for this run", type=str, required=True)
    parser.add_argument("--skip_alignment", action="store_true", help="skip alignment for this analysis, "
                                                                      "ONLY submit jobs for counting mutations in existing output folder")
    parser.add_argument("-r1", help="r1 SAM file", type=str)
    parser.add_argument("-r2", help="r2 SAM file", type=str)

    # user input arguments with default values set
    parser.add_argument("-log", "--log_level", help="set log level: debug, \
        info, warning, error, critical. (default = info)", type=str, default="info")
    parser.add_argument("-env", "--environment", help= "The cluster used to \
        run this script (default = DC)",type=str, default="DC")
    parser.add_argument("-at", type = int, help="Alignment time \
        (default = 8h)", default=8)
    parser.add_argument("-mt", type = int, help="Mutation call time \
        (default = 36h)", default=36)
    parser.add_argument("-mm", type=int, help="Mutation call request memory \
            (default = 15GB)", default=15)
    parser.add_argument("-c", type=int, help="Number of cores", default=16)
    parser.add_argument("-b", "--base", help="ASCII code base", default=33)
    parser.add_argument("-test", action="store_true", help="Turn on testing mode")
    parser.add_argument("-rc", action="store_true", help="Turn on rc mode, both direction of the reads will be "
                                                         "aligned to the reference. Variant calling will be "
                                                         "performed on all the reads that are aligned (BE CAREFUL!)")

    parser.add_argument("-override", "--sr_Override", action="store_true", help="Provide this argument when there is only one replicate")
    parser.add_argument("--posteriorQC", action="store_true", help="Turn on posterior QC mode, this requires more "
                                                                   "memory and runtime, please change the variable "
                                                                   "accordingly")
    parser.add_argument("--wt_override", action="store_true", help="When no wt conditions defined in the parameter sheet, turn on this option will treat EVERYTHING as wt. Phred scores will be adjusted based on the first replicate")

    parser.add_argument("--calibratePhredWT", action="store_true", help="When this parameter is provided, "
                                                                      "use wt to calibrate phred scores")
    parser.add_argument("--calibratePhredPhix", action="store_true", help="When this parameter is provided, "
                                                                      "use phix alignments to calibrate phred scores")
    parser.add_argument("--resubmit", action="store_true", help="For a finished run, batch resubmit failed scripts ("
                                                                "if any)")
    args = parser.parse_args()
    args.environment = args.environment.upper()

    # main(args)
