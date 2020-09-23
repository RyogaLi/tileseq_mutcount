#!/usr/bin/env python3.6

# Main script for sequencing analysis
# Input of this script:
#        - JSON file containing input parameters
#    - Fastq files (gzipped)

# what does this script do?
# 1. Read input csv file and convert to json (using csv2Json.R (you need to
# install required packages in R)
# 1. Read paramters from json file
# 2. Read input fastq file
#    2.a. Ramdonly select 30k reads and save a copy of downsampled fastq files

# 3. Align fastq file to reference sequence, generate sam files
# 4. From sam files, count mutations
# 5. Output mutation counts to summary.csv

# modules
import pandas as pd
import numpy as np
import os
import glob
import argparse
import logging
import datetime
import shutil
# import time

# pakage modules
import settings
import help_functions
import count_mut
import cluster


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
        self._main_path = os.path.abspath(__file__)  # path to this python file
        self._logging = main_log
        self._param_json = param_path  # parameter json file

        self._args = args  # user input arguments

        if args.fastq:
            self._fastq_list = glob.glob(args.fastq+"/*.fastq.gz")
        else:
            self._fastq_list = []
        self._output = output_dir

        # parse parameter json file
        self._project, self._seq, self._cds_seq, self._tile_map, self._region_map, \
            self._samples, self._var = help_functions.parse_json(param_path)
        self._sample_names = self._samples["Sample ID"].tolist()
        self._project = self._project.replace(" ", "_")  # project name
        # make main log
        self._log = self._logging.getLogger("main.log")


    def _init_skip(self, skip, r1=None, r2=None):
        """
        if skip == T: skip alignment
            make suboutput dir with time stamp
        else: not skip alignment
            make main output dir with time stamp
        """
        if skip:  # only run mutation counts on the aligned samples
            self._skip = True
            if args.r1 and args.r2:
                self._r1 = args.r1
                self._r2 = args.r2
                # self._log.info(f"SAM R1: {r1}")
                # self._log.info(f"SAM R2: {r2}")
            else:
                self._r1 = ""
                self._r2 = ""
                self._log.info(f"Sam files are read from {self._output}")
        else:
            self._r1 = ""
            self._r2 = ""
            self._skip = False

    def _align_sh_(self):
        """
        1. Align each fastq file to reference (submit one job for each alignment)
        2. Align each downsampled fastq file to reference (submit one job for each alignment)
        Output sam files into output folders (sam_files, ds_sam_files)
        """
        # VALIDATE: if all samples are present in fastq files
        # also check if any fastq file is empty
        fastq_sample_id = [os.path.basename(i).split("_")[0] for i in self._fastq_list]
        fastq_sample_id = list(set(fastq_sample_id))

        # validation
        # check if input fastq files contains all the samples in paramter file
        if set(self._sample_names).issubset(fastq_sample_id):
            self._log.info("All samples found")
            self._log.info(f"In total there are {len(list(set(self._sample_names)))} samples in the csv file")
            self._log.info(f"In total there are {len(fastq_sample_id)} fastq files")
        else:
            test= list(np.setdiff1d(self._sample_names,fastq_sample_id))
            join_list = ",".join(test)
            self._log.error("fastq files do not match input samples.")
            self._log.error(f"Fastq files not found for {join_list}")
            exit(1)

        # create mappving for R1 and R2 for each sample
        # ONLY samples provided in the parameter files will be analyzed
        fastq_map = []
        for r1 in self._fastq_list:
            ID = os.path.basename(r1).split("_")[0]
            if ("_R1_" in r1) and (ID in self._sample_names):
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
        if args.environment == "GURU":
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

        elif args.environment == "BC2" or args.environment == "DC" or args.environment == "BC":
            # make sh files to submit to BC
            sh_output = os.path.join(self._output, "BC_aln_sh")
            os.system("mkdir "+sh_output)

            if args.environment == "BC2" or args.environment == "BC":
                # make sh files to submit to BC
                self._log.info("Submitting alignment jobs to BC/BC2...")
                # alignment_sh_bc2(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging)
                sam_df, job_list = cluster.alignment_sh_bc2(fastq_map, self._project, self._seq.seq.values.item(), ref_path, sam_output, sh_output, args.at, self._log)
                self._log.info("Alignment jobs are submitte to BC2. Check pbs-output for STDOUT/STDERR")

            elif args.environment == "DC":
                # make sh files to submit to DC
                sh_output = os.path.join(self._output, "DC_aln_sh")
                os.system("mkdir "+sh_output)

                self._log.info("Submitting alignment jobs to DC...")
                sam_df, job_list = cluster.alignment_sh_dc(fastq_map, self._project, self._seq.seq.values.item(), ref_path, sam_output, sh_output, args.at, self._log)
                self._log.info("Alignment jobs are submitte to DC. Check pbs-output for STDOUT/STDERR")

            self._log.info(f"Total jobs submitted: {len(job_list)}")
            finished = cluster.parse_jobs(job_list, self._logging.getLogger("track.jobs"))  # track list of jobs

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
                            if i == 0:
                                line = line.split(" ")
                                f_info.append(line[0]) # number of reads in r1
                            elif i == 5:
                                line = line.split(" ")
                                f_info.append(line[0]) # overall alignment rate for r1
                            elif i == 6:
                                line = line.split(" ")
                                f_info.append(line[0]) # number of reads in r2
                            elif i == 11:
                                line = line.split(" ")
                                f_info.append(line[0]) # overall alignment rate for r2
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
        self._log.info("Counting mutations")
        # submit job with main.py -r1 and -r2
        # run main.py with -r1 and -r2

        mut_counts = count_mut.readSam(self._r1, self._r2, self._param_json, self._args, self._output)
        mut_counts._merged_main()

    def _makejobs(self, sh_output, sam_dir):
        """
        For each pair of sam files in output/sam_files/
        submit mut count job to the cluster
        """
        if sam_dir == "":
            # get samples in param file
            sam_dir = os.path.join(args.output, "sam_files/") # read sam file from sam_file
        # get sam files from parameter file
        if not os.path.isdir(sam_dir):
            self._log.error(f"Directory: ./sam_files/ not found in {self._output}")
            exit(1)
        self._log.debug(f"Sam files are read from {sam_dir}")
        job_list = []
        for i in self._sample_names:
            sam_f_r1 = glob.glob(f"{sam_dir}{i}_*R1_*.sam") # assume all the sam files have the same name format (id_*.sam)
            sam_f_r2 = glob.glob(f"{sam_dir}{i}_*R2_*.sam")
            if len(sam_f_r1) == 0 or len(sam_f_r2) == 0:
                self._log.error(f"SAM file for sample {i} not found. Please check your parameter file")
                exit(1)
            else:
                self._log.info(f"Sample {i}")
                self._log.info(f"Read1: {sam_f_r1[0]}")
                self._log.info(f"Read2: {sam_f_r2[0]}")
                sam_id_1 = os.path.basename(sam_f_r1[0]).split("_")[0]
                sam_id_2 = os.path.basename(sam_f_r2[0]).split("_")[0]
                if (sam_id_1 != i) or (sam_id_1 != sam_id_2) or (sam_id_2 != i):
                    self._log.error("IDs in sam files don't match!")
                self._r1 = sam_f_r1[0]
                self._r2 = sam_f_r2[0]

            # submit job with main.py -r1 and -r2
            # run main.py with -r1 and -r2
            if args.sr_Override:
                cmd = f"python {self._main_path} -n {args.name} -r1 {self._r1} -r2 {self._r2} -o {self._output} -p {self._param_json} --skip_alignment -log {args.log_level} -env {args.environment} -at {args.at} -mt {args.mt} -override"
            else:
                cmd = f"python {self._main_path} -n {args.name} -r1 {self._r1} -r2 {self._r2} -o {self._output} -p {self._param_json} --skip_alignment -log {args.log_level} -env {args.environment} -at {args.at} -mt {args.mt}"

            if args.environment == "BC2" or args.environment == "BC":
                logging.info("Submitting mutation counts jobs to BC2...")
                job_id = cluster.mut_count_sh_bc(i, cmd, args.mt, sh_output, self._log)
            elif args.environment == "DC":
                logging.info("Submitting mutation counts jobs to DC...")
                # (sample_name, cmd, mt, sh_output_dir, logger)
                job_id = cluster.mut_count_sh_dc(i, cmd, args.mt, sh_output, self._log) # this function will make a sh file for submitting the job

            job_list.append(job_id)

        jobs = ",".join(job_list)
        self._log.debug(f"All jobs: {jobs}")
        self._log.info(f"Total jobs running: {len(job_list)}")
        finished = cluster.parse_jobs(job_list, self._logging.getLogger("track.jobs")) # track list of jobs

        return finished

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


    def _main(self):
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
            if self._skip == False:
                time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
                # make mut count dir alignemnt not skipped
                # otherwise the user input dir should be the output dir gnerated
                self._output = os.path.join(self._output, args.name + "_" + time_now + "_mut_count")
                os.makedirs(self._output)

            # output directory is the mut_count dir
            # make folder to store all the sh files
            if args.environment == "BC2" or args.environment == "DC" or args.environment == "BC":
                sh_output = os.path.join(self._output, "BC_mut_sh")
                self._log.info(f"Mutation count sh files are made in {sh_output}")
                os.mkdir(sh_output)
            else:
                self._log.error("Please provide valid environment: BC/BC2/DC")
                exit(1)
            finished = self._makejobs(sh_output, sam_dir)

            if finished:
                self._log.info(f"Mutaion counting jobs are finished!")
                self._log.info("Check mutation counts file ...")
                mutcount_list = glob.glob(os.path.join(self._output, "counts_sample_*.csv"))
                self._log.info(f"{len(mutcount_list)} mutation counts file generated")
                if len(self._sample_names) > len(mutcount_list):
                    self._log.error("Job finished but some variant call files are missing!!")
                    missing = list(set(self._sample_names)-set(mutcount_list))
                    missing = ",".join(missing)
                    self._log.error(f"Missing samples: {missing}")
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
    # check if the correct args are provided
    # if you don't specify --skip_alignment then you cannot provide r1 and r2 sam_files
    if not args.skip_alignment:
        if args.r1 or args.r2:
            print(f"Invalid paramters! Please specify --skip_alignment if you want to analyze\
            one pair of sam files")
            exit(1)
    # check if output dir exists
    if not os.path.isdir(args.output):
        print(f"Output directory not found: {args.output}")
        exit(1)
    # check if fastq dir exists
    if args.fastq and not os.path.isdir(args.fastq):
        print(f"Fastq file path not found: {args.fastq}")
        exit(1)
    # try convert csv to json in the same dir as the csv file
    # convert csv file to json
    if args.param.endswith(".csv"):
        print("parameter sheet in csv format, converting to json")
        param_json = args.param.replace(".csv", ".json")
        if os.path.isfile(param_json):
             os.remove(param_json)
        if args.sr_Override:
            convert = f"Rscript {settings.CSV2JSON} {args.param} -o {param_json} --srOverride"
        else:
            convert = f"Rscript {settings.CSV2JSON} {args.param} -o {param_json}"
        os.system(convert)
    # if the file ends with .json, do nothing
    elif args.param.endswith(".json"):
        param_json = args.param
    else:
        print("Please provide valid paramter file format (csv or json)")
        exit(1)
    if not os.path.isfile(param_json):
        print("Json file does not exist, check conversion!")
        exit(1)

    return param_json

def main(args):
    """
    Main for fastq2counts
    """
    # get time stamp for this object
    time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    param_json = check(args)
    # get basename for the param file
    param_base = os.path.basename(param_json)
    if args.skip_alignment: # if we want to skip the alignment part
        # if user provided R1sam + R2sam
        if args.r1 and args.r2:
            # analyze one pair of r1 and r2 file
            # VALIDATE if both files are sam files
            if not args.r1.endswith(".sam") or not args.r2.endswith(".sam"):
                print("Please provide SAM files")
                return
            # copy the parameter file to output dir if the json file is not in the dir
            param_path = os.path.join(args.output, param_base)
            if not os.path.isfile(param_path):
                param_path = shutil.copy(param_json, args.output, follow_symlinks=True)
            # set up loggingthis creates main.log in the mut_count output dir
            main_log = log(args.output, args.log_level)

            # initialize mutcount object
            mc = fastq2counts(param_path, args.output, main_log, args)
            mc._init_skip(skip=True, r1=args.r1, r2=args.r2)

        else:
            # user provided a csv file and path to sam files
            # for each pair of sam file we would submit one job to the cluster for mutation counting
            # make time stamped output folder for this project
            updated_out = os.path.join(args.output, args.name + "_" + time_now + "_mut_count")
            # sam_dir = os.path.join(args.output, "sam_files/") # read sam file from sam_file
            os.makedirs(updated_out)
            args_log_path = os.path.join(updated_out, "args.log")
            write_param(args_log_path, args)
            # load json file
            param_path = os.path.join(updated_out, param_base)
            if not os.path.isfile(param_path):
                param_path = shutil.copy(param_json, updated_out, follow_symlinks=True)
            main_log = log(updated_out, args.log_level)
            # initialize mutcount object
            mc = fastq2counts(param_path, updated_out, main_log, args)
            mc._init_skip(skip=True)
    else:
        # alignment
        # make time stamped output folder for this project
        updated_out = os.path.join(args.output, args.name + "_" + time_now)
        os.makedirs(updated_out)  # make directory to save this run
        # write args to file
        args_log_path = os.path.join(updated_out, "args.log")
        write_param(args_log_path, args)

        param_path = os.path.join(updated_out, param_base)
        if not os.path.isfile(param_path):
            param_path = shutil.copy(param_json, updated_out, follow_symlinks=True)
        main_log = log(updated_out, args.log_level)

        # initialize mutcount object
        mc = fastq2counts(param_path, updated_out, main_log, args)
        mc._init_skip(skip=False)

    # start the run
    mc._main()

def write_param(args_log_path, args):
    """
    Write input arguments to param.log
    """
    with open(args_log_path, "w") as args_log:
        for arg in vars(args):
            args_log.write(arg+",")
            args_log.write(f"{getattr(args, arg, 'N/A')}\n")

def log(output_dir, log_level):
    """
    Make a logging object which writes to the main.log in output_dir
    """
    log_level = log_level.upper()
    # init main log
    logging.basicConfig(filename=os.path.join(output_dir, "main.log"),
                    filemode="a",
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                    datefmt="%m/%d/%Y %I:%M:%S %p",
                    level = log_level)
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(log_level)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s - %(name)-8s: %(levelname)-4s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    return logging


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='TileSeq mutation counts')
    # user input arguments
    parser.add_argument("-f", "--fastq", help="Path to all fastq files you want to analyze", type=str)
    parser.add_argument("-o", "--output", help="Output folder", type=str, required=True)
    parser.add_argument("-p", "--param", help="csv paramter file", type=str, required=True)
    parser.add_argument("--skip_alignment", action="store_true", help="skip alignment for this analysis, ONLY submit jobs for counting mutations in existing output folder")
    parser.add_argument("-n", "--name", help="Name for this run", type=str, required=True)
    parser.add_argument("-r1", help="r1 SAM file", type=str)
    parser.add_argument("-r2", help="r2 SAM file", type=str)

    # user input arguments with default values set
    parser.add_argument("-log", "--log_level", help="set log level: debug, \
        info, warning, error, critical. (default = debug)", type=str, default="debug")
    parser.add_argument("-env", "--environment", help= "The cluster used to \
        run this script (default = DC)",type=str, default="DC")
    ##parser.add_argument("-qual", "--quality", help="Posterior threshold for \
    ##    filtering mutations (default = 0.99)", type=float, default = 0.99)
    ##parser.add_argument("-min", "--min_cover", help="Minimal percentage required to \
    ##    cover the tile (default = 0.4)", type=float, default=0.6)
    parser.add_argument("-at", type = int, help="Alignment time \
        (default = 8h)", default=8)
    parser.add_argument("-mt", type = int, help="Mutation call time \
        (default = 36h)", default=36)
    parser.add_argument("-override", "--sr_Override", action="store_true", help="Provide this argument when there is only one replicate")

    args = parser.parse_args()

    print(" **** NOTE: Before you run this pipeline, please check settings.py to update program paths **** ")
    main(args)
