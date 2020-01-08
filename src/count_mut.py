#!~/usr/bin/python3

## Read sam file (R1 and R2)
# count mutations in sam files 
# output mutation counts

import pandas as pd
import numpy as np
import re
import os
import sys
import glob
import pysam
import logging
import argparse
import datetime

from pathlib import Path

# modules in package 
import help_functions
import locate_mut

class readSam(object):

    def __init__(self, sam_r1, sam_r2, seq_lookup, tile_map, region_map, samples, output_dir, qual_filter, log_level):
        """
        sam_R1: read one of the sample
        sam_R2: read two of the sample
        seq_lookup: df contains DNA sequence, Protein sequence and positions mapped to each other

        tile_map: tile start and end pos
        region_map: region start and end pos

        samples: df with sample names and corresponding conditions, tiles
        output_dir: main output directory

        qual_filter: quality fileter number to sam file

        log level: settings for logging
    
        """
        self._sample_id = os.path.basename(sam_r1).split("_")[0]

        self._r1 = sam_r1
        self._r2 = sam_r2

        self._seq_lookup = seq_lookup

        self._tile_map = tile_map
        self._region_map = region_map
        
        self._qual = qual_filter

        time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        # create time stamped dir for mutation count output (if not created before)
        self._output_counts_dir = output_dir + time + "_mut_call"
        try:
            os.makedirs(self._output_counts_dir)
        except FileExistsError:
            pass

        # create logging file (time stamped)
        log_dir = os.path.join(self._output_counts_dir, "mut_log")
        try:
            os.makedirs(log_dir)
        except FileExistsError:
            # directory already exists
            pass

        # create time stamped log file for this mutation count call in the mut_log folder
        log_f = "sample_"+ str(self._sample_id)+"_"+"mut_count.log"
        
        logging.basicConfig(filename=os.path.join(log_dir, log_f), filemode="w", format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p", level = log_level.upper())
        logging.info(f"Counting mutations in sample-{self._sample_id}")
        logging.info(f"Sam file input R1: {sam_r1}")
        logging.info(f"Sam file input R2: {sam_r2}")


    def _convert_sam(self, input_sam):
        """
        trim sam file, convert it to a df with the following colmns:
        read_name, mapped_name, pos_start, mapQ, CIGAR, seq
        
        input sam files are filtered at this stage:
        * mapeed name should match input fasta sequence name 
        * mapQ > mapQ cut off value
        * mdz indicates mutations in the read

        """
        read_counts = 0
        trimmed = []
        # read sam file
        with open(input_sam, "r") as sam_file:
            for line in sam_file:

                if line.startswith("@"): # skip header line
                    continue
                
                read_counts += 1
                line = line.split("\t")
                # create table to save information from sam file
                read_name = line[0]

                mapped_name = line[2]
                if mapped_name == "*": # if the read didn't mep 
                    continue

                pos_start = line[3]

                mapQ = int(line[4])
                if mapQ < self._qual: # remove reads with mapQ < quality score filter
                    continue

                CIGAR = line[5]
                seq = line[9]
                quality = line[10]
                
                mdz = [i for i in line if "MD:Z:" in i]
                if len(mdz) != 0:
                    mdz = mdz[0].split(":")[-1]
                else:
                    mdz = ""
                if (not re.search('[a-zA-Z]', mdz)) and ("I" not in CIGAR): 
                    # remove reads that have no mutations in MDZ
                    continue
                trimmed.append([read_name, mapped_name, pos_start, mapQ, CIGAR, mdz, seq, quality])
        
        trimmed = pd.DataFrame(trimmed, columns=["read_name", "mapped_name", "pos_start", "mapQ", "CIGAR", "mdz","seq", "qual"])
        logging.info(f"Raw sequencing depth for {input_sam} (before filtering): {read_counts}")
        return trimmed

    def _main(self, full_seq, cds_seq):
        """
        """
        r1_df = self._convert_sam(self._r1)
        r2_df = self._convert_sam(self._r2)
        logging.info(f"Sequencing depth for R1: {r1_df.shape[0]}")
        logging.info(f"Sequencing depth for R2: {r2_df.shape[0]}")
        # merge r1 and r2 based on read ID
        merge = pd.merge(r1_df, r2_df, how="left", on="read_name", suffixes=("_r1", "_r2"))
        merge = merge.dropna()
        logging.info(f"After merging & filtering sam R1 and sam R2, {merge.shape[0]} read-pairs remained for analysis")
        # facts['pop2050'] = facts.apply(lambda row: final_pop(row['population'],row['population_growth']),axis=1)
        mut_reads=0 
        for index, row in merge.iterrows():
            # for each read pair in sam file, identify the mutation
            # process R1 and R2 together
            mut_parser = locate_mut.MutParser(row, full_seq, cds_seq, self._seq_lookup) 
            mp_update = mut_parser._get_seq()
            #mp_update = mut_parser._build_lookup()
            
            # get mutation from R1
            r1_mut = mp_update._parse_mut(mp_update._r1_cigar, mp_update._r1_mdz, mp_update._r1_ref, mp_update._r1_read, mp_update._r1_pos, mp_update._r1_qual)
            
            # get mutation from R2
            r2_mut = mp_update._parse_mut(mp_update._r2_cigar, mp_update._r2_mdz, mp_update._r2_ref, mp_update._r2_read, mp_update._r2_pos, mp_update._r2_qual)

            # get overlap of mutations in R1 and mutations in R2
            both_mut = list(set(r1_mut) & set(r2_mut))
            if not both_mut == []:
                mut_reads +=1
                mp_update._translate_mut(both_mut)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='TileSeq mutation counts (for sam files)') 
    parser.add_argument("-r1", "--read_1", help="sam file for R1", required=True)
    parser.add_argument("-r2", "--read_2", help="sam file for R2", required=True)
    parser.add_argument("-qual", "--quality", help="sam file mapQ filter", default=20)
    parser.add_argument("-o", "--output", help="Output folder", required = True)
    parser.add_argument("-log", "--log_level", help="set log level: debug, info, warning, error, critical.", default = "debug") 
    parser.add_argument("-p", "--param", help="json paramter file", required = True)
    parser.add_argument("-env", "--environment", help= "The cluster used to run this script", default="GURU")
    args = parser.parse_args()
    
    sam_r1 = args.read_1
    sam_r2 = args.read_2
    qual_filter = args.quality

    out = args.output
    param = args.param
    env = args.environment

    # process the json file 
    project, seq, cds_seq, tile_map, region_map, samples = help_functions.parse_json(param)
    # build lookup table
    lookup_df = help_functions.build_lookup(seq.cds_start.item(), seq.cds_end.item(), cds_seq)
    
    # initialize the object
    MutCounts = readSam(sam_r1, sam_r2, lookup_df, tile_map, region_map, samples, out, qual_filter, args.log_level)

    MutCounts._main(seq, cds_seq)

