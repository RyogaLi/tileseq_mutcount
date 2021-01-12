#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')

import os
import pandas as pd

# for each pair of sam files, check how many of the reads are flipped
d = sys.argv[1]
output = []
for f in os.listdir(d):

    if ("_R1_" in f) and (".sam" in f):

        fname = f.split("_R1_")[0]
        r2 = f"{fname}_R2_001.sam"

        r1_f = open(os.path.join(d, f), "r")
        r2_f = open(os.path.join(d, r2), "r")
        # init objects
        total_reads = 0
        flipped = 0
        r1_flipped = 0
        r2_flipped = 0
        for line_r1, line_r2 in zip(r1_f, r2_f):
            total_reads +=1
            line_r1 = line_r1.split()
            line_r2 = line_r2.split()
            # check if read ID mapped
            read_name_r1 = line_r1[0]
            read_name_r2 = line_r2[0]
            if read_name_r1 != read_name_r2:
                print(f"Read pair IDs did not map, please check fastq files {fname}")
                exit(1)
            # check flag, for R1, check for 16, for R2, check for 0
            if (int(line_r1[1]) == 16) and (int(line_r2[1]) == 0):
                # this pair is flipped
                flipped +=1
            if int(line_r1[1]) == 16:
                r1_flipped +=1

            if int(line_r2[1]) == 0:
                r2_flipped +=1
        output.append([fname, total_reads, flipped, r1_flipped, r2_flipped])
        print(output)
        print(fname)

df = pd.DataFrame(output, columns=["sample", "total reads", "flipped pairs", "r1 flipped", "r2 flipped"])
print(df)
df.to_csv("./count_rc.csv", index=False)



