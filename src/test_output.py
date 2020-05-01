## load a output file and check hgvs

import pandas as pd
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
from hgvs.easy import parser, projector

def check(mutcount, nm):
    """
    mutcount = file generated from the mutation counts
    nm = NM id of the target gene
    """
    df = pd.read_csv(mutcount, skiprows=18)
    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    for index, row in df.iterrows():
        var_c = f"{nm}:{row['HGVS']}"
        var_c = hp.parse_hgvs_variant(var_c)
        print(var_c)
        var_p = hgvs.easy.am38.c_to_p(var_c)
        print(var_p)

if __name__ == "__main__":

    # MTHFR
    mutcount = "/home/home7/rothlab/rli/dev/tileseq_dev/output/MTHFR_new_alignm_2020-04-22-21-18-44/mut_count_test/counts_sample_47.csv"
    nm = "NM_005957.4"

    # SUMO

    check(mutcount, nm)
