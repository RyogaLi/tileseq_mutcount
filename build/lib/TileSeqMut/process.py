#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys

sys.path.append('..')

import os
import pandas as pd

import multiprocessing as mp

def process_wrapper(combined_dict):
  """
  Call another module
  """
  return combined_dict

def main(f1, f2, cores):

  # open files and read lines
  file1 = open(f1, "r")
  file2 = open(f2, "r")

  # init pool objects
  pool = mp.Pool(cores)
  jobs = []
  combined = {}

  for line_f1, line_f2 in zip(file1, file2):
    # some preprocess of the lines
    line_f1 = line_f1.split()
    line_f2 = line_f2.split()

    if line_f1[0] != line_f2[0]:
      exit(1)

    if "#" in line_f1:
      continue

    # add some information from both lines into a dictionary
    combined["ID_f1"] = line_f1[0]
    combined["ID_f2"] = line_f2[0]

    combined["qual_f1"] = line_f1[3]
    combined["qual_f2"] = line_f2[3]

    # pass this dictionary to the process
    jobs.append(pool.apply_async(process_wrapper, (combined)))

  file1.close()
  file2.close()

  # get output from all the jobs and store the output
  id_list, name_list = [], []
  for job in jobs:
    ids, names = job.get()
    id_list.append(ids)
    name_list.append(names)

  # clean up
  pool.close()

if __name__ == '__main__':
  f1 = "/home/rothlab2/rli/tileseq_output/201011_LDLR_Surface_2020-10-15-17-54-25/sam_files/LDLR-Wildtype-T06-Rep2_90k_R1_001.sam"
  f2 = "/home/rothlab2/rli/tileseq_output/201011_LDLR_Surface_2020-10-15-17-54-25/sam_files/LDLR-Wildtype-T06-Rep2_90k_R2_001.sam"

  main(f1, f2, 8)