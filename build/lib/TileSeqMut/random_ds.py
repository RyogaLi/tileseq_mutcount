#!/usr/bin/python3.7
# to run this script
# python random_sample.py -input /path/to/fastq.gz
# - all the files in the input directory have to be gzipped
# - the script output down sampled fastq files in the same directory as the original fastq files
# to move them to a new dir you can run `mv *_ds.fastq /new_dir/`
# output fastq files are NOT gzipped, you will need to zip them inorder to run song's script

import pandas as pd
import os
import random
import argparse


def random_lines(r1, r2, n):
    """
    get n random records from input file
    output to file *_ds.fastq
    """
    record_number = 0

    base_r1 = os.path.basename(r1)
    base_r2 = os.path.basename(r2)

    output_r1 = os.path.dirname(r1)+"/"+base_r1.replace(".fastq", "_ds.fastq")
    output_r2 = os.path.dirname(r2)+"/"+base_r2.replace(".fastq", "_ds.fastq")

    with open(r1, "r") as i:
        r1_lines = sum([1 for line in i])
    with open(r2, "r") as j:
        r2_lines = sum([1 for line in j])

    total_records_r1 = int(r1_lines/4)
    total_records_r2 = int(r2_lines/4)

    if n > total_records_r1 or n > total_records_r2:
        os.system("cp "+r1+ " "+output_r1)
        os.system("cp "+r2+ " "+output_r2)
        return None
    else:
        if total_records_r1 <= total_records_r2:
            records_to_keep = set(random.sample(range(total_records_r1 + 1), n))
        else:
            records_to_keep = set(random.sample(range(total_records_r2 + 1), n))

    with open(r1, "r") as i_1:
        with open(output_r1, "w") as o1:
            for line1 in i_1:
                line2 = i_1.readline()
                line3 = i_1.readline()
                line4 = i_1.readline()
                if record_number in records_to_keep:
                    o1.write(line1)
                    o1.write(line2)
                    o1.write(line3)
                    o1.write(line4)
                record_number += 1

    record_number = 0
    with open(r2, "r") as i_2:
        with open(output_r2, "w") as o2:
            for line1 in i_2:
                line2 = i_2.readline()
                line3 = i_2.readline()
                line4 = i_2.readline()
                if record_number in records_to_keep:
                    o2.write(line1)
                    o2.write(line2)
                    o2.write(line3)
                    o2.write(line4)
                record_number += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-input", help="input folder contains original fastq files")
    parser.add_argument("-n", help="number of lines to downsample to")

    args = parser.parse_args()

    n = args.n

    for f in os.listdir(args.input):
        if f.endswith(".fastq.gz") and ("_R1_" in f):

            r2 = f.replace("_R1_", "_R2_")

            r1 = os.path.join(args.input, f)
            r2 = os.path.join(args.input, r2)

            os.system("gunzip "+r1)
            os.system("gunzip "+r2)

            r1 = r1.replace(".gz", "")
            r2 = r2.replace(".gz", "")

            random_lines(r1, r2, n)
