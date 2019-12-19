#!/usr/bin/python3.6

# Main script for sequencing analysis
# Input of this script: 
#   - JSON file containing input parameters
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
import argparse
import glob
import sys
import json
import logging
import datetime

# pakage modules
import settings

class MutCount(object):

    def __init__(self, json_file, output_folder):
        """
        Initialize mutation counts 
        Load input json file and parameters for this run
        json_file: path to param.json
        output_folder: user input folder for storing time stamped output folders
        From param.json load all the parameters, save them into df
        """
        
        # check if output folder exists
        if not os.path.isdir(out):
            print(f"Output directory not found: {out}")
            exit(1)
        
        # get time stamp for this object
        time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

        with open(json_file, "r") as jsonf:
            data = json.load(jsonf)
            self._project = data["project"] # project name
            
            # make time stamped output folder for this project
            self._output = self._project.replace(" ", "-")
            self._output = os.path.join(output_folder, self._output+ "_" +time)
            
            mkdir_cmd = "mkdir "+self._output
            print(mkdir_cmd)

            # create main log file in this output folder
            logging.basicConfig(filename=os.path.join(self._output, "main.log"))
            
            self._tmp_seq = data["template"]["seq"] # template sequence
            self._cds_seq = self._tmp_seq # cds sequence
            
            # validation:
            if len(self._cds_seq) % 3 != 0: 
                pass
            
            self._tile_map = ""
            self._region_map = ""



def submit(jobs):
    """
    function to submit jobs to sge
    jobs: list of sh files with command
    """
    for j in jobs:
        basename = os.path.basename(j)
        job_name = basename.split(".")[0]
        cmd = "qsub -cwd -V -N "+job_name + " "+j
        os.system(cmd) 

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='TileSeq mutation counts') 
    parser.add_argument("-f", "--fastq", help="Path to all fastq files you want to analyze")
    parser.add_argument("-o", "--output", help="Output folder", required = True)
    parser.add_argument("-d", "--debug", help="printing to the terminal instead of log file") 
    parser.add_argument("-p", "--params", help="path to paramter file", required = True)

    args = parser.parse_args()

    f = args.fastq
    out = args.output
    param = args.params

    # Initialize MutCount main 
    mc = MutCount(param, out)

    
    
