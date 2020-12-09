#!/usr/bin/env python3.6
# testing script for all the python scripts
import pandas as pd
import os
import argparse
import timeit
import logging

import count_mut as cm

# test different scripts in the package

class UnitTest(object):

    """Unitesting for all the python scripts"""

    def __init__(self, args):
        self.args = args

    def _test_count_mut(self):
        """
        Test count_mut.py
        """
        # init all the input parameters
        # sam_r1, sam_r2, param, log_level, output_dir
        sam_r1 = "/home/rothlab1/rli/dev/tilseq_mutcount/output/SUMO_rerun_2020-05-01-15-59-22/sam_files/30_R1_joint.sam"
        sam_r2 = "/home/rothlab1/rli/dev/tilseq_mutcount/output/SUMO_rerun_2020-05-01-15-59-22/sam_files/30_R2_joint.sam"
        param = "/home/rothlab1/rli/dev/tilseq_mutcount/output/SUMO_rerun_2020-05-01-15-59-22/SUMO1_parameters.json"
        log_level = "DEBUG"
        output_dir = "/home/rothlab1/rli/dev/tilseq_mutcount/output/SUMO_rerun_2020-05-01-15-59-22/testdel_2020-07-14-12-43-26_mut_count"
        # init read sam
        print("Testing count_mut")
        start = timeit.default_timer()
        test_readSam = cm.readSam(sam_r1, sam_r2, param, log_level, output_dir, 1)
        test_readSam.multi_core()
        stop = timeit.default_timer()
        # without any processing, simply go through the files as chunks
        # Linear time to read both sam files: 10.06786861596629
        # without mutation calling, but with other constructs
        # Linear time to read both sam files: 169.94309764588252
        print(f"Linear time to read both sam files: {stop - start}")

    def _test_loc_mut(self):
        """
        """
        pass

def main(args):
    ut = UnitTest(args)
    if args.m == "count_mut":
        ut._test_count_mut()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Unit Test modules')
    # user input arguments
    parser.add_argument("-m", help="Mode of testing")
    args = parser.parse_args()
    main(args)
