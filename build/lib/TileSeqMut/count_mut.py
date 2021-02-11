#!/usr/bin/env python3.7

## Read sam file (R1 and R2)
# count mutations in sam files
# output mutation counts

import pandas as pd
import re
import os
import math
import logging
import mmap
import multiprocessing as mp
import argparse
# import datetime
from collections import deque

# modules in package
from TileSeqMut import help_functions
from TileSeqMut import locate_mut

# import help_functions
# import locate_mut

class readSam(object):

    def __init__(self, sam_r1, sam_r2, param, arguments, output_dir, cores, loggerobj):
        """
        sam_R1: read one of the sample
        sam_R2: read two of the sample

        output_dir: main output directory

        log_level: settings for logging
        log_dir: directory to save all the logging for each sample
        """
        self._r1 = sam_r1
        self._r2 = sam_r2
        self._project, self._seq, self._cds_seq, self._tile_map, self._region_map, self._samples, self._var = help_functions.parse_json(param)

        self._qual = self._var["posteriorThreshold"]
        min_cover = self._var["minCover"]
        self._mutrate = self._var["mutRate"]
        self._output_counts_dir = output_dir
        self._cores = cores

        # sample information
        self._sample_id = os.path.basename(sam_r1).split("_")[0]
        self._sample_info = self._samples[self._samples["Sample ID"] == self._sample_id]

        # tile information
        self._sample_tile = self._sample_info["Tile ID"].values[0]
        self._tile_begins = (self._tile_map[self._tile_map["Tile Number"] == self._sample_tile]["Start AA"].values[0] *3)-2 # beginning position of this tile (cds position)
        self._tile_ends = self._tile_map[self._tile_map["Tile Number"] == self._sample_tile]["End AA"].values[0] *3  # ending position of this tile (cds position)
        self._tile_len = self._tile_ends - self._tile_begins
        self._cds_start = self._seq.cds_start
        self._min_map_len = math.ceil(self._tile_len * float(min_cover))

        self._sample_condition = self._sample_info["Condition"].values[0]
        self._sample_tp = self._sample_info["Time point"].values[0]
        self._sample_rep = self._sample_info["Replicate"].values[0]

        self._seq_lookup = help_functions.build_lookup(self._cds_start.values.item(), self._seq.cds_end.values.item(), self._cds_seq)

        # init a new object for logging, log to sample specific file.log
        # log_f = os.path.join(output_dir, f"sample_{str(self._sample_id)}_mut_count.log")
        # log_object = help_functions.logginginit(arguments.log_level, log_f)

        self._sample_counts_f = os.path.join(self._output_counts_dir,f"counts_sample_{self._sample_id}.csv")
        self._sample_counts_r1_f = os.path.join(self._output_counts_dir, f"counts_sample_{self._sample_id}_r1.csv")
        self._sample_counts_r2_f = os.path.join(self._output_counts_dir, f"counts_sample_{self._sample_id}_r2.csv")

        self._base = arguments.base
        # self._mut_log = log_object.getLogger("count.mut")
        # self._locate_log = log_object.getLogger("locate.mut")
        self._mut_log = loggerobj
        self._mut_log.info(f"Counting mutations in sample-{self._sample_id}")
        self._mut_log.info(f"Sam file input R1:{sam_r1}")
        self._mut_log.info(f"Sam file input R2:{sam_r2}")

        output_csv = open(self._sample_counts_f, "w")
        # write log information to counts output
        output_csv.write(f"#Sample:{self._sample_id}\n#Tile:{self._sample_tile}\n#Tile Starts:{self._tile_begins}\n"
                         f"#Tile Ends:{self._tile_ends}\n#Condition:{self._sample_condition}\n"
                         f"#Replicate:{self._sample_rep}\n#Timepoint:{self._sample_tp}\n"
                         f"#Posterior cutoff:{self._qual}\n#min cover %:{min_cover}\n")
        output_csv.close()

        output_csv = open(self._sample_counts_r1_f, "w")
        # write log information to counts output
        output_csv.write(f"#Sample:{self._sample_id}\n#Tile:{self._sample_tile}\n#Tile Starts:{self._tile_begins}\n"
                         f"#Tile Ends:{self._tile_ends}\n#Condition:{self._sample_condition}\n"
                         f"#Replicate:{self._sample_rep}\n#Timepoint:{self._sample_tp}\n"
                         f"#Posterior cutoff:{self._qual}\n#min cover %:{min_cover}\n")
        output_csv.close()

        output_csv = open(self._sample_counts_r2_f, "w")
        # write log information to counts output
        output_csv.write(f"#Sample:{self._sample_id}\n#Tile:{self._sample_tile}\n#Tile Starts:{self._tile_begins}\n"
                         f"#Tile Ends:{self._tile_ends}\n#Condition:{self._sample_condition}\n"
                         f"#Replicate:{self._sample_rep}\n#Timepoint:{self._sample_tp}\n"
                         f"#Posterior cutoff:{self._qual}\n#min cover %:{min_cover}\n")
        output_csv.close()


    def _merged_main(self):
        """
        Read two sam files at the same time, store mutations that passed filter
        """
        read_pair = 0 # total pairs
        un_map = 0 # total number of unmapped reads
        read_nomut = 0 # read pairs that have no mutations

        hgvs_output = {} # save all the hgvs in this pair of sam files
        off_mut = {} # save all the positions that were mapped but outside of the tile
        row = {} # create a dictionary to save all the reads information to pass on to locate_mut.py

        final_pairs = 0
        off_read = 0

        chunkSize = 1500000 # number of characters in each chunk (you will need to adjust this)
        chunk1 = deque([""]) #buffered lines from 1st file
        chunk2 = deque([""]) #buffered lines from 2nd file
        with open(self._r1, "r") as r1_f, open(self._r2, "r") as r2_f:
            while chunk1 and chunk2:
                line_r1 = chunk1.popleft()
                if not chunk1:
                    line_r1,*more = (line_r1+r1_f.read(chunkSize)).split("\n")
                    chunk1.extend(more)
                line_r2 = chunk2.popleft()
                if not chunk2:
                    line_r2,*more = (line_r2+r2_f.read(chunkSize)).split("\n")
                    chunk2.extend(more)

                if line_r1.startswith("@") or line_r2.startswith("@"): # skip header lines
                    continue

                line_r1 = line_r1.split()
                line_r2 = line_r2.split()
                if len(line_r1) <= 1:
                    continue

                read_pair += 1 # count how many read pairs in this pair of sam files
                mapped_name_r1 = line_r1[2]
                mapped_name_r2 = line_r2[2]
                if mapped_name_r1 == "*" or mapped_name_r2 == "*": # if one of the read didn't map to ref
                    un_map +=1
                    continue

                # check if read ID mapped
                read_name_r1 = line_r1[0]
                read_name_r2 = line_r2[0]
                if read_name_r1 != read_name_r2:
                    print(read_name_r1)
                    print(read_name_r2)
                    self._mut_log.error("Read pair IDs did not map, please check fastq files")
                    exit(1)

                # get starting position for r1 and r2
                pos_start_r1 = line_r1[3]
                pos_start_r2 = line_r2[3]
                # record reads that are not mapped to this tile
                # this is defined as if the read covers at least min_map percent of the tile
                # the default min_map is 70%
                # we also assume that the read len is ALWAYS 150
                r1_end = int(pos_start_r1) + 150
                # r1_end must cover from start of the tile to 70% of the tile
                if (r1_end - int(self._cds_start)) < (int(self._tile_begins) + int(self._min_map_len)):
                    off_read += 1
                    continue
                if (int(pos_start_r2) - int(self._cds_start)) > (int(self._tile_ends) - int(self._min_map_len)):
                    off_read += 1
                    continue
                # get CIGAR string
                CIGAR_r1 = line_r1[5]
                seq_r1 = line_r1[9]
                quality_r1 = line_r1[10]

                CIGAR_r2 = line_r2[5]
                seq_r2 = line_r2[9]
                quality_r2 = line_r2[10]

                mdz_r1 = [i for i in line_r1 if "MD:Z:" in i]
                mdz_r2 = [i for i in line_r2 if "MD:Z:" in i]

                if len(mdz_r1) != 0:
                    mdz_r1 = mdz_r1[0].split(":")[-1]
                else:
                    mdz_r1 = ""

                if len(mdz_r2) != 0:
                    mdz_r2 = mdz_r2[0].split(":")[-1]
                else:
                    mdz_r2 = ""

                if ((not re.search('[a-zA-Z]', mdz_r1)) and ("I" not in CIGAR_r1)) and ((not re.search('[a-zA-Z]', mdz_r2)) and ("I" not in CIGAR_r2)):
                    # if MDZ string only contains numbers
                    # and no insertions shown in CIGAR string
                    # means there is no mutation in this read
                    # if both reads have no mutations in them, skip this pair
                    read_nomut +=1
                    # remove reads that have no mutations in MDZ
                    continue

                # make the reads in the format of a dictionary
                # columns=["mapped_name", "pos_start", "qual", "CIGAR", "mdz","seq"])

                #row["mapped_name_r1"] = mapped_name_r1
                row["pos_start_r1"] = pos_start_r1
                row["qual_r1"] = quality_r1
                row["cigar_r1"] = CIGAR_r1
                row["mdz_r1"] = mdz_r1
                row["seq_r1"] = seq_r1

                #row["mapped_name_r2"] = mapped_name_r2
                row["pos_start_r2"] = pos_start_r2
                row["qual_r2"] = quality_r2
                row["cigar_r2"] = CIGAR_r2
                row["mdz_r2"] = mdz_r2
                row["seq_r2"] = seq_r2

                # pass this dictionary to locate mut
                # mut = locate_mut_main()
                # add mutation to mut list
                mut_parser = locate_mut.MutParser(row, self._seq, self._cds_seq, self._seq_lookup,
                                                  self._tile_begins, self._tile_ends, self._qual,
                                                  self._mut_log, self._mutrate, self._base)
                hgvs, outside_mut, pos_df, all_posterior= mut_parser._main()
                if len(hgvs) !=0:
                    final_pairs +=1
                    if hgvs in hgvs_output:
                        hgvs_output[hgvs] += 1
                    else:
                        hgvs_output[hgvs] = 1
                    #hgvs_output.append(hgvs)
                if outside_mut != []:
                    outside_mut = list(set(outside_mut))
                    for i in outside_mut:
                        if not (i in off_mut):
                            off_mut[i] = 1

        # track mutations that are not within the same tile
        # track sequencing depth for each sample
        # write this information to a tmp file
        # merge the tmp file (in main.py)
        tmp_f = os.path.join(self._output_counts_dir,f"{self._sample_id}_tmp.csv")
        off_mut = list(set(off_mut))
        off_mut.sort()
        f = open(tmp_f, "w")
        f.write(f"sample,{self._sample_id},tile,{self._tile_begins}-{self._tile_ends},sequencing_depth,{read_pair},off_tile_reads,{off_read},off_tile_perc,{off_read/read_pair}\n")
        f.close()

        output_csv = open(self._sample_counts_f, "a")
        output_csv.write(f"#Raw read depth:{read_pair}\n")
        output_csv.write(f"#Number of read pairs without mutations:{read_nomut}\n#Number of read pairs did not map to gene:{un_map}\n#Mutation positions outside of the tile:{off_mut}\n#Number of reads outside of the tile:{off_read}\n")
        output_csv.write(f"#Final read-depth:{read_pair - un_map - off_read}\n")
        output_csv.write(f"#Total read pairs with mutations:{final_pairs}\n")
        output_csv.write(f"#Comment: Total read pairs with mutations = Read pairs with mutations that passed the posterior threshold\n#Comment: Final read-depth = raw read depth - reads didn't map to gene - reads mapped outside of the tile\n")

        output_csv.close()

        self._mut_log.info(f"Raw sequencing depth: {read_pair}")
        self._mut_log.info(f"Number of reads without mutations:{read_nomut}")
        self._mut_log.info(f"Final read-depth: {read_pair - un_map - off_read}")
        # convert list to df with one col
        hgvs_df = pd.DataFrame.from_dict(hgvs_output, orient="index")
        hgvs_df = hgvs_df.reset_index()
        if not hgvs_df.empty:
            hgvs_df.columns = ["HGVS", "count"]
        else:
            hgvs_df = pd.DataFrame({}, columns=["HGVS", "count"])
        hgvs_df.to_csv(self._sample_counts_f, mode="a", index=False)

    def multi_core(self):
        """
        Read two sam files at the same time, store mutations that passed filter
        """
        read_pair = 0  # total pairs
        un_map = 0  # total number of unmapped reads
        read_nomut = 0  # read pairs that have no mutations

        hgvs_output = {}  # save all the hgvs in this pair of sam files
        r1_pop_hgvs = {}  # save all the hgvs that are 2/3nt changes based only on R1
        r2_pop_hgvs = {}  # save all the hgvs that are 2/3nt changes based only on R2
        off_mut = {}  # save all the positions that were mapped but outside of the tile
        row = {}  # create a dictionary to save all the reads information to pass on to locate_mut.py

        final_pairs = 0
        r1_popmut = 0
        r2_popmut = 0
        off_read = 0

        # chunkSize = 1500000  # number of characters in each chunk
        # chunk1 = deque([""])  # buffered lines from 1st file
        # chunk2 = deque([""])  # buffered lines from 2nd file
        r1_f = open(self._r1, "r")
        r2_f = open(self._r2, "r")
        # init objects
        pool = mp.Pool(self._cores)
        jobs = []
        self._mut_log.info("Start reading SAM files")
        for line_r1, line_r2 in zip(r1_f, r2_f):

            line_r1 = line_r1.split()
            line_r2 = line_r2.split()
            if len(line_r1) < 9 or len(line_r2) < 9:
                # the read has no sequence
                self._mut_log.warning("Missing fields in read!")
                continue

            read_pair += 1  # count how many read pairs in this pair of sam files
            mapped_name_r1 = line_r1[2]
            mapped_name_r2 = line_r2[2]
            if mapped_name_r1 == "*" or mapped_name_r2 == "*":  # if one of the read didn't map to ref
                un_map += 1
                continue

            # check if read ID mapped
            read_name_r1 = line_r1[0]
            read_name_r2 = line_r2[0]
            if read_name_r1 != read_name_r2:
                self._mut_log.error("Read pair IDs did not map, please check fastq files")
                exit(1)

            # get starting position for r1 and r2
            pos_start_r1 = line_r1[3]
            pos_start_r2 = line_r2[3]
            # record reads that are not mapped to this tile
            # this is defined as if the read covers at least min_map percent of the tile
            # the default min_map is 70%
            # we also assume that the read len is ALWAYS 150
            r1_end = int(pos_start_r1) + 150
            # r1_end must cover from start of the tile to 70% of the tile
            if (r1_end - int(self._cds_start)) < (int(self._tile_begins) + int(self._min_map_len)):
                off_read += 1
                continue
            if (int(pos_start_r2) - int(self._cds_start)) > (int(self._tile_ends) - int(self._min_map_len)):
                off_read += 1
                continue

            # get CIGAR string
            CIGAR_r1 = line_r1[5]
            seq_r1 = line_r1[9]
            quality_r1 = line_r1[10]

            CIGAR_r2 = line_r2[5]
            seq_r2 = line_r2[9]
            quality_r2 = line_r2[10]

            mdz_r1 = [i for i in line_r1 if "MD:Z:" in i]
            mdz_r2 = [i for i in line_r2 if "MD:Z:" in i]

            if len(mdz_r1) != 0:
                mdz_r1 = mdz_r1[0].split(":")[-1]
            else:
                mdz_r1 = ""

            if len(mdz_r2) != 0:
                mdz_r2 = mdz_r2[0].split(":")[-1]
            else:
                mdz_r2 = ""

            if (not re.search('[a-zA-Z]', mdz_r1)) and ("I" not in CIGAR_r1) and ("D" not in CIGAR_r2) and (not re.search('[a-zA-Z]', mdz_r2)) and ("I" not in CIGAR_r2) and ("D" not in CIGAR_r2):
                # if MDZ string only contains numbers
                # and no insertions shown in CIGAR string
                # means there is no mutation in this read
                # if both reads have no mutations in them, skip this pair
                read_nomut += 1
                # remove reads that have no mutations in MDZ
                continue
            # make the reads in the format of a dictionary
            # columns=["mapped_name", "pos_start", "qual", "CIGAR", "mdz","seq"])
            # row["mapped_name_r1"] = mapped_name_r1
            row["pos_start_r1"] = pos_start_r1
            row["qual_r1"] = quality_r1
            row["cigar_r1"] = CIGAR_r1
            row["mdz_r1"] = mdz_r1
            row["seq_r1"] = seq_r1

            # row["mapped_name_r2"] = mapped_name_r2
            row["pos_start_r2"] = pos_start_r2
            row["qual_r2"] = quality_r2
            row["cigar_r2"] = CIGAR_r2
            row["mdz_r2"] = mdz_r2
            row["seq_r2"] = seq_r2
            # mut_parser = locate_mut.MutParser(row, self._seq, self._cds_seq, self._seq_lookup, self._tile_begins,
                                              # self._tile_ends, self._qual, self._locate_log, self._mutrate)
            jobs.append(pool.apply_async(process_wrapper, (row, self._seq, self._cds_seq, self._seq_lookup, self._tile_begins,
                                               self._tile_ends, self._qual, self._mut_log, self._mutrate,
                                                           self._base)))
            row = {} # flush

        r1_f.close()
        r2_f.close()
        self._mut_log.info("File streamed to subprocesses, waiting for jobs to finish")
        # wait for all jobs to finish
        all_df = []

        for job in jobs:
            hgvs, outside_mut, all_posterior, hgvs_r1_clusters, hgvs_r2_clusters = job.get()
            if len(hgvs) != 0:
                final_pairs += 1
                if hgvs in hgvs_output:
                    hgvs_output[hgvs] += 1
                else:
                    hgvs_output[hgvs] = 1
                # hgvs_output.append(hgvs)

            if len(hgvs_r1_clusters) != 0:
                r1_popmut += 1
                if hgvs_r1_clusters in r1_pop_hgvs:
                    r1_pop_hgvs[hgvs_r1_clusters] += 1
                else:
                    r1_pop_hgvs[hgvs_r1_clusters] = 1

            if len(hgvs_r2_clusters) != 0:
                r2_popmut += 1
                if hgvs_r2_clusters in r2_pop_hgvs:
                    r2_pop_hgvs[hgvs_r2_clusters] += 1
                else:
                    r2_pop_hgvs[hgvs_r2_clusters] = 1

            if outside_mut != []:
                outside_mut = list(set(outside_mut))
                for i in outside_mut:
                    if not (i in off_mut):
                        off_mut[i] = 1
            all_df.append(all_posterior)

        # clean up
        pool.close()
        self._mut_log.info("Pool closed")
        # track mutations that are not within the same tile
        # track sequencing depth for each sample
        # write this information to a tmp file
        # merge the tmp file (in main.py)

        # tmp_f = os.path.join(self._output_counts_dir, f"{self._sample_id}_tmp.csv")
        # off_mut = list(set(off_mut))
        # off_mut.sort()
        # f = open(tmp_f, "w")
        # f.write(
        #     f"sample,{self._sample_id},tile,{self._tile_begins}-{self._tile_ends},sequencing_depth,{read_pair},off_tile_reads,{off_read},off_tile_perc,{off_read / read_pair}\n")
        # f.close()

        output_csv = open(self._sample_counts_f, "a")
        output_csv.write(f"#Raw read depth:{read_pair}\n")
        output_csv.write(
            f"#Number of read pairs without mutations:{read_nomut}\n#Number of read pairs did not map to gene:{un_map}\n#Mutation positions outside of the tile:{off_mut}\n#Number of reads outside of the tile:{off_read}\n")
        output_csv.write(f"#Final read-depth:{read_pair - un_map - off_read}\n")
        output_csv.write(f"#Total read pairs with mutations:{final_pairs}\n")
        output_csv.write(
            f"#Comment: Total read pairs with mutations = Read pairs with mutations that passed the posterior threshold\n#Comment: Final read-depth = raw read depth - reads didn't map to gene - reads mapped outside of the tile\n")

        output_csv.close()

        self._mut_log.info(f"Raw sequencing depth: {read_pair}")
        self._mut_log.info(f"Number of reads without mutations:{read_nomut}")
        self._mut_log.info(f"Final read-depth: {read_pair - un_map - off_read}")
        # convert list to df with one col
        hgvs_df = pd.DataFrame.from_dict(hgvs_output, orient="index")
        hgvs_df = hgvs_df.reset_index()
        if not hgvs_df.empty:
            hgvs_df.columns = ["HGVS", "count"]
        hgvs_df.to_csv(self._sample_counts_f, mode="a", index=False)

        r1_df = pd.DataFrame.from_dict(r1_pop_hgvs, orient="index")
        r1_df = r1_df.reset_index()
        if not r1_df.empty:
            r1_df.columns = ["HGVS", "count"]
        r1_df.to_csv(self._sample_counts_r1_f, mode="a", index=False)

        r2_df = pd.DataFrame.from_dict(r2_pop_hgvs, orient="index")
        r2_df = r2_df.reset_index()
        if not r2_df.empty:
            r2_df.columns = ["HGVS", "count"]
        r2_df.to_csv(self._sample_counts_r2_f, mode="a", index=False)
        #
        # self._mut_log.info(f"Reading posterior filtered dfs ...")
        #
        # # pos_df = pd.concat(pos_df)
        # posterior_f = os.path.join(self._output_counts_dir, f"{self._sample_id}_posprob.csv")
        # for df in pos_df:
        #     df.to_csv(posterior_f, mode='a', header=False, index=False)

        self._mut_log.info(f"Reading posterior unfiltered dfs ...")
        # all_df = pd.concat(all_df)
        all_f = os.path.join(self._output_counts_dir, f"{self._sample_id}_posprob_all.csv")
        for df in all_df:
            df.to_csv(all_f, mode='a', header=False, index=False)
        self._mut_log.info(f"Posterior df saved to files ... {all_f}")


def process_wrapper(row, seq, cds_seq, seq_lookup, tile_begins, tile_ends, qual, locate_log, mutrate, base):
    """

    """
    mut_parser = locate_mut.MutParser(row, seq, cds_seq, seq_lookup, tile_begins, tile_ends, qual, locate_log,
                                      mutrate, base)
    hgvs, outside_mut, all_df, hgvs_r1_clusters, hgvs_r2_clusters = mut_parser._main()
    return hgvs, outside_mut, all_df, hgvs_r1_clusters, hgvs_r2_clusters


# if __name__ == "__main__":
#
#     parser = argparse.ArgumentParser(description='TileSeq mutation counts (for sam files)')
#     parser.add_argument("-r1", "--read_1", help="sam file for R1", required=True)
#     parser.add_argument("-r2", "--read_2", help="sam file for R2", required=True)
#     parser.add_argument("-qual", "--quality", help="quality filter using posterior probability", default=0.99)
#     parser.add_argument("-o", "--output", help="Output folder that contains the directory: ./sam_files", required = True)
#     parser.add_argument("-log", "--log_level", help="Set log level: debug, info, warning, error, critical.", default = "debug")
#     parser.add_argument("-mutlog", "--log_dir", help="Directory to save all the log files for each sample")
#     parser.add_argument("-p", "--param", help="csv paramter file", required = True)
#     parser.add_argument("-min", "--min_cover", help="Minimum percentage required to be covered by reads", default = 0.4)
#
#     args = parser.parse_args()
#
#     sam_r1 = args.read_1
#     sam_r2 = args.read_2
#     qual_filter = float(args.quality)
#     log_dir = args.log_dir
#     min_map = float(args.min_cover)
#
#     out = args.output
#     param = args.param
#
#     # conver the csv file to json
#     # csv2json = os.path.abspath("src/csv2json.R")
#     param_json = param.replace(".csv", ".json")
#     #convert = f"Rscript {csv2json} {param} -o {param_json} -l stdout"
#     #os.system(convert)
#
#     # process the json file
#     project, seq, cds_seq, tile_map, region_map, samples = help_functions.parse_json(param_json)
#     # build lookup table
#     lookup_df = help_functions.build_lookup(seq.cds_start.values.item(), seq.cds_end.values.item(), cds_seq)
#
#     # initialize the object
#     MutCounts = readSam(sam_r1, sam_r2, seq, lookup_df, tile_map, region_map, samples, out, qual_filter, min_map, args.log_level, log_dir)
#
#     MutCounts._merged_main(seq, cds_seq)
