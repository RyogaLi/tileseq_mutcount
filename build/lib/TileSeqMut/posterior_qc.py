#!/usr/bin/env python3.7

import sys
import os
import pandas as pd
from ast import literal_eval
import seaborn as sns
import matplotlib.pyplot as plt
sys.path.append('..')

from TileSeqMut import help_functions
# import help_functions

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

# QC script for posteriors

class PosteriorQC(object):

    def __init__(self, mut_count_dir, param_json):
        """
        @param mut_count_dir: Directory ends with _mut_count
        @param param_json: Json parameter sheet
        """
        self._project, self._seq, self._cds_seq, self._tile_map, \
        self._region_map, self._samples, self._var = help_functions.parse_json(param_json)

        self._mut_dir = mut_count_dir

    def plot_counts(self, df):
        """
        @param df: all_df or pos_df
        """
        pass

    def runQC(self):
        """
        Run posterior QC, generate posterior QC plots in a subdir
        """
        # list all the prob files in the _mut_count dir
        for f in os.listdir(self._mut_dir):
            full_f_path = os.path.join(self._mut_dir, f)
            if "posprob.csv" in f:
                full_f_path = os.path.join(self._mut_dir, f)
                # get sample ID
                sample_id = f.split("_")[0]
                # match sample ID to tile
                sample_info = self._samples[self._samples["Sample ID"] == sample_id]
                sample_tile = sample_info["Tile ID"].values[0]

                prob_filtered = pd.read_csv(full_f_path)
                full_all_path = os.path.join(self._mut_dir, f"{sample_id}_posprob_all.csv")
                prob_all = pd.read_csv(full_all_path)



def read_pos_file(pos_file, all_file):
    """
    @param pos_file: file contains all the mutations and calculated posteriors for each mutation

    """
    pos_df = pd.read_csv(all_file)

    # mutations only on read 1
    read1_only = pos_df.loc[pos_df["read"] == "r1"]
    # mutations only on read 2
    read2_only = pos_df.loc[pos_df["read"] == "r2"]

    # mutations found on both reads
    both_read = pos_df.loc[(pos_df["read"] != "r1") & (pos_df["read"] != "r2")]
    both_read['read'] = [literal_eval(x) for x in both_read['read']]
    both_read['prob'] = [literal_eval(x) for x in both_read['prob']]
    both_read[["prob_r1", "prob_r2"]] = pd.DataFrame(both_read['prob'].tolist(), index=both_read.index)
    both_read[["label_r1", "label_r2"]] = pd.DataFrame(both_read['read'].tolist(), index=both_read.index)

    # make plot for mutation on r1 only and r2 only vs all
    n_r1 = read1_only.shape[0]
    n_r2 = read2_only.shape[0]
    both = both_read.shape[0]

    bar_df = pd.DataFrame({"type": ["Mutations on read1", "Mutations on read2", "Mutations on both reads"],
                           "count": [n_r1, n_r2, both]})
    sns.barplot(x="type", y="count", data=bar_df, palette=["#F07D78", "#B7F86E", "#316FA3"])
    plt.title("Number of mutations (before filtering) on each read")
    plt.show()
    plt.close()

    pos_df = pd.read_csv(pos_file)
    # mutations only on read 1
    read1_only = pos_df.loc[pos_df["read"] == "r1"]
    # mutations only on read 2
    read2_only = pos_df.loc[pos_df["read"] == "r2"]

    # mutations found on both reads
    both_read = pos_df.loc[(pos_df["read"] != "r1") & (pos_df["read"] != "r2")]
    if not both_read.empty:
        both_read['read'] = [literal_eval(x) for x in both_read['read']]
        both_read['prob'] = [literal_eval(x) for x in both_read['prob']]
        both_read[["prob_r1", "prob_r2"]] = pd.DataFrame(both_read['prob'].tolist(), index=both_read.index)
        both_read[["label_r1", "label_r2"]] = pd.DataFrame(both_read['read'].tolist(), index=both_read.index)
        both = both_read.shape[0]
    else:
        both = 0
    n_r1 = read1_only.shape[0]
    n_r2 = read2_only.shape[0]

    bar_df = pd.DataFrame({"type": ["Mutations on read1", "Mutations on read2", "Mutations on both reads"],
                           "count": [n_r1, n_r2, both]})
    sns.barplot(x="type", y="count", data=bar_df, palette=["#F07D78", "#B7F86E", "#316FA3"])
    plt.title("Number of mutations (before filtering) on each read")
    plt.show()
    plt.close()

if __name__ == '__main__':
    all_file = "/Users/roujia/Desktop/200819_WAS_QC_all_2020-08-19-13-11-35/56_posprob_all.csv"
    pos_file = "/Users/roujia/Desktop/200819_WAS_QC_all_2020-08-19-13-11-35/56_posprob.csv"

    # all_file = "/Users/roujia/Desktop/200819_WAS_QC_all_2020-08-19-13-11-35/56_posprob_all.csv"
    # pos_file = "/Users/roujia/Desktop/200819_WAS_QC_all_2020-08-19-13-11-35/56_posprob.csv"
    read_pos_file(pos_file, all_file)