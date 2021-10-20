#!/usr/bin/env python3.7

# Helper functions
# 1. downsample fastq files into n reads
# 2. loading process animation function
# 3. parse input json file

import pandas as pd
import random
import os
import sys
import json
import time
import itertools
import logging
from Bio.Seq import Seq


def parse_json(json_file):
    """
    parse input json file
    return information about the run
    """
    with open(json_file, "r") as jsonf:
        data = json.load(jsonf)
        project = data["project"] # project name
        seq = pd.DataFrame(data["template"], index=[0]) # sequence
        cds_seq = seq.seq[0]
        cds_seq = cds_seq[int(seq.cds_start)-1: int(seq.cds_end)]
        # create dictionary to map tile/region # to start,end positions
        try:
            tile_map = pd.DataFrame(data["tiles"])
        except:
            tile_map = pd.DataFrame(data["tiles"], index=[0])

        try:
            region_map = pd.DataFrame(data["regions"])
        except:
            region_map = pd.DataFrame(data["regions"], index=[0])
        try:
            samples = pd.DataFrame.from_dict(data["samples"])
        except:
            samples = pd.DataFrame.from_dict(data["samples"], index=[0])

        relations = data["conditions"]["definitions"]
        # get posterior quality cut off
        var_caller = data["varcaller"]
    return project, seq, cds_seq, tile_map, region_map, samples, var_caller, relations


def loadingAnimation(process):
    """
    loading animation to show when loading a process
    """

    while process.is_alive():
        chars = "/—\|"
        for char in chars:
            sys.stdout.write('\r'+'downsamping reads ... '+char)
            time.sleep(.1)
            sys.stdout.flush()
    sys.stdout.write("\n")


def build_lookup(cds_start, cds_end, cds_seq):
    """
    build a lookup df for a given gene
    return a df with columns:
    temp_pos, cds, dna_pos, protein, protein_pos
    """
    lookup_table = {}
    # list of template pos
    temp = list(range(cds_start, cds_end+1))
    # list of coding DNA bases
    cDNA = [i for i in cds_seq]
    cDNA_pos = list(range(1, len(cDNA)+1))
    # list of protein bases
    protein = [i for i in Seq(cds_seq).translate()]
    protein_pos = range(1, len(protein)+1)
    lookup_table["temp_pos"] = temp
    lookup_table["cds"] = cDNA
    lookup_table["cds_pos"] = cDNA_pos
    lookup_table["protein"] = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in protein))
    lookup_table["protein_pos"] = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in protein_pos))
    df = pd.DataFrame(lookup_table)
    return df


def logginginit(log_level, log_f):
    """
    Init logging in console as well as main log
    """
    log_level = log_level.upper()
    logging.basicConfig(filename=log_f,
                        filemode="w",
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                        datefmt="%m/%d/%Y %I:%M:%S %p",
                        level=log_level)

    # define a Handler which writes INFO messages or higher to the sys.stdout
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(log_level)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s - %(name)-8s: %(levelname)-4s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # stderr_logger = logging.getLogger('STDERR')
    # sl = StreamToLogger(stderr_logger, logging.ERROR)
    # sys.stderr = sl

    return logging


class StreamToLogger(object):
   """
   Fake file-like stream object that redirects writes to a logger instance.
   """
   def __init__(self, logger, log_level=logging.INFO):
      self.logger = logger
      self.log_level = log_level
      self.linebuf = ''

   def write(self, buf):
      for line in buf.rstrip().splitlines():
         self.logger.log(self.log_level, line.rstrip())


if __name__ == "__main__":
    j = "/home/rothlab2/rli/tileseq_output/210416_CHEK2_QC_2021-04-19-16-02-03/210416_optimizedRef_CHEK2.json"
    parse_json(j)
