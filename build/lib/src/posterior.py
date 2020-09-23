#!/usr/bin/env python3.6

# This is used to calculate the posterior probability for variants
# the script is modified from varCallerSnippet.r

import math
from fractions import Fraction

def parse_mut_cluster(mut_cluster, mut_rate):
    """
    Parse cluster of mutations
    mut_cluster is a  list of clusters (clustered mutations found in one read pair)
    Mutations are clustered together if they are within
    """
    pass

def bayesian_variant_call(basecall, phred, wt, mut_rate, clusterSize=1):
    """
    basecall: list of base calls (i.e R1 -> A R2 -> C :  ["A", "C"])
    phred: phred score for the base calls (in letters) ["!", "J"]
    wt: wild type base
    mut_rate: mutation rate

    return: dictinary with basecall as keys and post prob as values
    """
    nt = ["A", "G", "C", "T"] # all possible hypo bases
    # convert phred to int scores
    phred = [10**(-(ord(i) - 33) / 10) for i in phred]
    post_p = []
    for base in nt: # go through each nt
        log_odd = 0
        if base == wt:
            log_odd += math.log((1-mut_rate) ** clusterSize) - math.log(1-((1-mut_rate)**clusterSize)/3)
        else:
            log_odd += math.log(1-(1-mut_rate) ** clusterSize) - math.log(3) - math.log(1-(1-(1-mut_rate) ** clusterSize)/3)

        # insertion prior
        # log_odd += math.log(mut_rate) - math.log(4) - math.log(1-(mut_rate/4))

        for j in range(len(basecall)):
            if basecall[j] == base:
                log_odd += (math.log(1-phred[j]) - math.log(phred[j]) + math.log(3))
            else:
                log_odd += (math.log(phred[j]) - math.log(3) - math.log((1/3) -(phred[j]/9)))

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
