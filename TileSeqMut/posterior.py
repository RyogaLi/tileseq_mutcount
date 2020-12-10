#!/usr/bin/env python3.7

# This is used to calculate the posterior probability for variants
# the script is modified from varCallerSnippet.r

import math
import pandas as pd
import numpy as np
from fractions import Fraction

def cluster(mut_cluster:dict, mut_rate:float, cut_off:float):
    """
    Parse cluster of mutations
    mut_cluster is a  list of clusters (clustered mutations found in one read pair)
    Mutations are clustered together if they are within
    @param mut_cluster dictionary where the values are dfs with mutations from both r1 and r2
    @param mut_rate mut rate defined by user
    @param cut_off posterior cutoff defined by user
    """
    pos_prob = {"m": [], "prob": [], "read": []}
    for c in mut_cluster.keys():
        mutcall = mut_cluster[c]
        # each mutcall is a df with columns:
        # m_r1,read,pos,ref_r1,alt_r1,qual_r1,m_r2,ref_r2,alt_r2,qual_r2
        c_size = mutcall.shape[0]  # size of the cluster
        # iterate through the clusters

        for index, row in mutcall.iterrows():

            if (pd.isnull(row["m_r1"]) or row["alt_r1"] == "N") and not pd.isnull(row["ref_r2"]):
                pos = bayesian_variant_call([row["alt_r2"]], [row["qual_r2"]], row["ref_r2"], mut_rate, c_size)

                if pos[next(iter(pos))] > cut_off:
                    pos_prob["m"].append(row["m_r2"])
                    pos_prob["prob"].append(list(pos.values())[0])
                    pos_prob["read"].append("r2")

            elif (pd.isnull(row["m_r2"]) or row["alt_r2"] == "N") and not pd.isnull(row["ref_r1"]):
                pos = bayesian_variant_call([row["alt_r1"]], [row["qual_r1"]], row["ref_r1"], mut_rate, c_size)
                if pos[next(iter(pos))] > cut_off:
                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append(list(pos.values())[0])
                    pos_prob["read"].append("r1")

            elif (not pd.isnull(row["m_r2"])) and (not pd.isnull(row["ref_r1"])):


                basecall = [row["alt_r1"], row["alt_r2"]]
                qual = [row["qual_r1"], row["qual_r2"]]
                pos = bayesian_variant_call(basecall, qual, row["ref_r1"], mut_rate, c_size)
                if pos[row["alt_r1"]] > pos[row["alt_r2"]] and pos[row["alt_r1"]] > cut_off:

                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append(pos[row["alt_r1"]])
                    pos_prob["read"].append("r1")
                elif pos[row["alt_r2"]] > pos[row["alt_r1"]] and pos[row["alt_r2"]] > cut_off:
                    pos_prob["m"].append(row["m_r2"])
                    pos_prob["prob"].append(pos[row["alt_r2"]])
                    pos_prob["read"].append("r2")

                elif pos[row["alt_r1"]] == pos[row["alt_r2"]] and pos[row["alt_r1"]] > cut_off:
                    pos_prob["m"].append(row["m_r2"])
                    pos_prob["prob"].append(pos[row["alt_r2"]])
                    pos_prob["read"].append("-")

    # print(pos_prob)
    pos_df = pd.DataFrame(pos_prob)

    return pos_df


def bayesian_variant_call(basecall, qual, wt, mut_rate, clusterSize=1):
    """
    basecall: list of base calls (i.e R1 -> A R2 -> C :  ["A", "C"])
    phred: phred score for the base calls (in letters) ["!", "J"]
    wt: wild type base
    mut_rate: mutation rate

    return: dictinary with basecall as keys and post prob as values
    """
    # all possible hypo bases

    nt = list(set([wt]+basecall))
    # nt = ["A", "G", "C", "T"]
    # convert phred to int scores
    # convert phred to int scores
    # in the case of deletions or insertions, there are multiple phred scores for each mut call
    phred = []
    for i in qual:
        if len(i) == 1:
            phred.append(10 ** (-(ord(i) - 33) / 10))
        else:
            all_phred = [10 ** (-(ord(j) - 33) / 10) for j in i.split(",")]
            phred.append(np.prod(all_phred))

    # phred = [10**(-(ord(i) - 33) / 10) for i in phred]
    post_p = []
    for base in nt: # go through each nt
        log_odd = 0
        if base == wt:
            log_odd += math.log((1-mut_rate) ** clusterSize) - math.log(1-((1-mut_rate)**clusterSize)/3)
        elif base == "ins":
            log_odd += (math.log(mut_rate) - math.log(4)) - math.log(1-(mut_rate/4))
        else:
            log_odd += math.log(1-(1-mut_rate) ** clusterSize) - math.log(3) - math.log(1-(1-(1-mut_rate) ** clusterSize)/3)

        # insertion prior
        # log_odd += math.log(mut_rate) - math.log(4) - math.log(1-(mut_rate/4))

        for j in range(len(basecall)):
            if basecall[j] == base:
                log_odd += (math.log(1-phred[j]) - math.log(phred[j]) + math.log(3))
            else:
                log_odd += (math.log(phred[j]) - math.log(3) - math.log((1/3) -(phred[j]/9)))
        # print(log_odd)
        # print(basecall)
        logit_value = math.exp(log_odd) / (1+math.exp(log_odd))
        post_p.append(logit_value)

    prob = dict(zip(nt, post_p))

    output = dict(zip(basecall, [prob.get(base) for base in basecall]))

    return output

if __name__ == "__main__":
    basecall = ["T", "A"]
    phred = ["I", "J"]
    wt = "C"
    mut_rate= 0.0025
    prob = bayesian_variant_call(basecall, phred, wt, mut_rate)
    print(prob)
