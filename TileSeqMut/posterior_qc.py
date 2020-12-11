#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')
import os
import pandas as pd
from ast import literal_eval
import seaborn as sns
import matplotlib.pyplot as plt
# QC script for posteriors

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
    all_file = "/Users/roujia/Desktop/200819_WAS_QC_all_2020-08-19-13-11-35/80_posprob_all.csv"
    pos_file = "/Users/roujia/Desktop/200819_WAS_QC_all_2020-08-19-13-11-35/80_posprob.csv"
    read_pos_file(pos_file, all_file)