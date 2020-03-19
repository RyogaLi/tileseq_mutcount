#~/lib/Python-3.6.4/python

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
import os
import glob
import argparse
import logging
import datetime
import shutil
import time

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

        self._fastq_list = glob.glob(args.fastq+"/*.fastq.gz")
        self._output = output_dir

        # parse parameter json file
        self._project, self._seq, self._cds_seq, self._tile_map, self._region_map, \
            self._samples = help_functions.parse_json(param_path)
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
                self._log.info(f"SAM R1: {r1}")
                self._log.info(f"SAM R2: {r2}")
            else:
                self._r1 = ""
                self._r2 = ""
                self._log.info(f"Sam files are read from {self._output}")
        else:
            self._skip = False

    def _align_sh_(self):
        """
        1. Align each fastq file to reference (submit one job for each alignment)
        2. Align each downsampled fastq file to reference (submit one job for each alignment)
        Output sam files into output folders (sam_files, ds_sam_files)
        """
        # VALIDATE: if all samples are present in fastq files
        # also check if any fastq file is empty

        sample_names = self._samples["Sample ID"].tolist()
        fastq_sample_id = [os.path.basename(i).split("_")[0] for i in self._fastq_list]
        fastq_sample_id = list(set(fastq_sample_id))

        # validation
        # check if input fastq files contains all the samples in paramter file
        if set(sample_names).issubset(fastq_sample_id):
            self._log.info("All samples found")
            self._log.info(f"In total there are {len(list(set(sample_names)))} samples in the csv file")
            self._log.info(f"In total there are {len(fastq_sample_id)} fastq files")
        else:
            self._log.error("fastq files do not match input samples.")
            self._log.error("Program terminated due to error")
            exit(1)

        # create mappving for R1 and R2 for each sample
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

        mut_counts = count_mut.readSam(self._r1, self._r2, self._param, self._args, self._output)
        mut_counts._merged_main()

    def _makejobs(self, sh_output):
        """
        For each pair of sam files in output/sam_files/
        submit mut count job to the cluster
        """
        # get samples in parameter file
        sample_names = self._samples["Sample ID"].tolist() # samples in param file
        # get sam files from parameter file
        sam_dir = os.path.join(self._output, "sam_files/") # read sam file from sam_file
        for i in sample_names:
            sam_f_r1 = glob.glob(f"{sam_dir}{i}_R1*.sam") # assume all the sam files have the same name format (id_*.sam)
            sam_f_r2 = glob.glob(f"{sam_dir}{i}_R2*.sam")
            if len(sam_f_r1) == 0 or len(sam_f_r2) == 0:
                logging.error(f"SAM file for sample {i} not found. Please check your parameter file")
                exit(1)
            else:
                self._r1 = sam_f_r1
                self._r2 = sam_f_r2

            # submit job with main.py -r1 and -r2
            # run main.py with -r1 and -r2
            cmd = f"python {self._main_path} \
                    -r1 {self._r1} -r2 {self._r2} -o {self._output} -p {self._param} --skip_alignment \
                    -log {args._log_level} -env {args.environment} -qual {args.qual} -min {args._min_cover} \
                    -at {args.at} -mt {args.mt}"
            job_list = []
            if args.environment == "BC2" or args.environment == "BC":
                logging.info("Submitting mutation counts jobs to BC2...")
                job_id = cluster.mut_count_sh_bc(i, cmd, sh_output, self._mt, self._log)
            elif args.environment == "DC":
                logging.info("Submitting mutation counts jobs to DC...")
                job_id = cluster.mut_count_sh_dc(i, cmd, sh_output, self._mt, self._log) # this function will make a sh file for submitting the job
            job_list.append(job_id)
        self._log.info(f"Total jobs running: {len(job_list)}")
        finished = cluster.parse_jobs(job_list, self._logging.getLogger("track.jobs")) # track list of jobs

        if finished:
            self._log.info(f"Mutaion counting jobs are finished!")
            self._log.info("Check mutation counts file ...")
            mutcount_list = glob.glob(os.path.join(self._output, "counts_sample_*.csv"))
            self._log.info(f"{len(mutcount_list)} mutation counts file generated")
            for f in mutcount_list:
                mut_n = 0
                # double check if each output file has mutations
                with open(f, "r") as mut_output:
                    for line in mut_output:
                        # skip header
                        if "c." in line:
                            mut_n += 1
                if mut_n == 0:
                    self._log.error(f"{f} has 0 variants! Check mut log for this sample.")
                else:
                    self._log.info(f"{f} has {mut_n} variants")

    def _main(self):
        """
        """
        if self._skip == False: # submit jobs for alignment
            # get sam_df from self._align_sh()
            self._align_sh_()

        if self._r1 != "" and self._r2 != "":
            # call functions in count_mut.py
            # init(sam_r1, sam_r2, param, args, ouptut_dir, logger)
            self._mut_count()
        # submit jobs for mutation counting
        # if user did not provide r1 and r2 SAM file
        # get r1 and r2 sam files from output dir provided by user
        elif self._r1 == "" and self._r2 == "":
            if not args.skip_alignment:
                time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
                # make mut count dir alignemnt not skipped
                # otherwise the user input dir should be the output dir gnerated
                self._output = os.path.join(args.output, args.name + "_" + time_now + "_mut_count")
                os.makedirs(self._output)

            # output directory is the mut_count dir
            # make folder to store all the sh files
            if args.environment == "BC2" or args.environment == "DC" or args.environment == "BC":
                self._log.info(f"Running on {self._env}")
                sh_output = os.path.join(self._output, "BC_mut_sh")
                os.mkdir(sh_output)
            # make folder to store all the log files
            log_dir = os.path.join(self._output, "mut_log")
            os.makedirs(log_dir)
            # get samples in parameter file

            self._makejobs(self, sh_output, self._output)

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
    if not os.path.isdir(args.fastq):
        print(f"Fastq file path not found: {args.fastq}")
        exit(1)
    # try convert csv to json in the same dir as the csv file
    # convert csv file to json
    if args.param.endswith(".csv"):
        param_json = args.param.replace(".csv", ".json")
        convert = f"Rscript {settings.CSV2JSON} {args.param} -o {param_json} --srOverride"
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

def main(args):
    """
    Main for fastq2counts
    Validate user inputs
    """
    # get time stamp for this object
    time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    param_json = check(args)
    # get basename for the param file
    param_base = os.path.basename(param_json)
    output_dir = ""
    if args.skip_alignment: # if we want to skip the alignment part
        # if user provided R1sam + R2sam
        if args.r1 and args.r2:
            # analyze one pair of r1 and r2 file
            # VALIDATE if both files are sam files
            if not args.r1.endswith(".sam") or not args.r2.endswith(".sam"):
                print("Please provide SAM files")
                exit(1)
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
            os.makedirs(updated_out)
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
        param_path = os.path.join(updated_out, param_base)
        if not os.path.isfile(param_path):
            param_path = shutil.copy(param_json, updated_out, follow_symlinks=True)
        main_log = log(updated_out, args.log_level)

        # initialize mutcount object
        mc = fastq2counts(param_path, updated_out, main_log, args)
        mc._init_skip(skip=False)

    if output_dir != "":
        # this parameter file will be saved in the main output dir
        param_f = open(os.path.join(output_dir, "param.log"), "a")
        # log - time for this run (same as the output folder name)
        param_f.write(f"Run started: {time}\n")
        param_f.write(f"Run name: {args.name}\n")
        param_f.write(f"This run was on the cluster: {args.env}\n")
        param_f.write(f"Input parameter file used: {args.param}\n")
        param_f.write(f"Posterior probability cutoff for mutation counts: {args.qual}\n")
        param_f.write(f"Min_ cover percentage: {args.min_cover}\n")
        param_f.close()

    # start the run
    mc._main()

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
    formatter = logging.Formatter('%(name)-8s: %(levelname)-4s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    return logging


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='TileSeq mutation counts')
    # user input arguments
    parser.add_argument("-f", "--fastq", help="Path to all fastq files you want to analyze")
    parser.add_argument("-o", "--output", help="Output folder", required=True)
    parser.add_argument("-p", "--param", help="csv paramter file", required=True)
    parser.add_argument("--skip_alignment", action="store_true", help="skip alignment for this analysis, ONLY submit jobs for counting mutations in existing output folder")
    parser.add_argument("-r1", help="r1 SAM file")
    parser.add_argument("-r2", help="r2 SAM file")

    # user input arguments with default values set
    parser.add_argument("-log", "--log_level", help="set log level: debug, \
        info, warning, error, critical. (default = debug)", default="debug")
    parser.add_argument("-env", "--environment", help= "The cluster used to \
        run this script (default = DC)", default="DC")
    parser.add_argument("-qual", "--quality", help="Posterior threshold for \
        filtering mutations (default = 0.99)", default = 0.99)
    parser.add_argument("-name", help="Name for this run", required=True)
    parser.add_argument("-min", "--min_cover", help="Minimal % required to \
        cover the tile (default = 0.4)", default=0.6)
    parser.add_argument("-at", type = int, help="Alignment time \
        (default = 8h)", default=8)
    parser.add_argument("-mt", type = int, help="Mutation call time \
        (default = 36h)", default=36)

    args = parser.parse_args()

    # convert all the args to variables
    # store all the parameters in dictionary
    # param_dict["f"] = args.fastq # path to fastq files
    # param_dict["out"] = args.output # path to ouptut dir
    # param_dict["param"] = args.param # path to parameter.csv
    # param_dict["env"] = args.environment # environment settings (DC, BC, BC2)
    # param_qual = float(args.quality) # posterior quality cutoff
    # min_cover = float(ais.min_cover) # min coverage
    # log_level = args.log_level # log level defaultl
    # r1sam = args.r1 # sam file r1
    # r2sam = args.r2 # sam file r2
    # mt = args.mt # mutation count time required
    # at = args.at # alignment time required
    # name = args.name
    print(" **** NOTE: Before you run this pipeline, please check settings.py to update program paths **** ")
    main(args)
