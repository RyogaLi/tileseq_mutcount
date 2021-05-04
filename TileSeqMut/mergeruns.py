#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import os
import pandas as pd
import json
import argparse
import subprocess

from TileSeqMut.help_functions import logginginit

#for two runs merge the mutation counts files into output folder

def merge_runs(input_dir1, input_dir2, map_df, output):
    """
    :param input_dir1: Directory 1 contains mutation counts files
    :param input_dir2: Directory 2 contains mutation counts files
    :param map_df: A dataframe that maps the samples to be merged
    :param output: output directory to save merged mut_counts
    """
    # make a new file

    new_samples = []
    for index, row in map_df.iterrows():
        merged_header = ""
        mut_count_f1 = os.path.join(input_dir1, f"counts_sample_{row['Sample ID_run1']}.csv")
        mut_count_f2 = os.path.join(input_dir2, f"counts_sample_{row['Sample ID_run2']}.csv")

        new_samples.append(f"{row['Sample ID_run1']}-{row['Sample ID_run2']}")
        # join headers
        with open(mut_count_f1) as file1, open(mut_count_f2) as file2:
            for line1, line2 in zip(file1, file2):
                if "#" not in line1: break
                if "Final read-depth:" in line1:
                    n_depth_1 = int(line1.split(":")[1])
                    n_depth_2 = int(line2.split(":")[1])
                    total_depth = n_depth_1+n_depth_2
                    merged_header+= f"#Final read-depth:{total_depth}\n"
                    continue
                if "Number of read pairs without mutations:" in line1:
                    n_depth_1 = int(line1.split(":")[1])
                    n_depth_2 = int(line2.split(":")[1])
                    total_depth = n_depth_1 + n_depth_2
                    merged_header += f"#Number of read pairs without mutations:{total_depth}\n"
                    continue
                if "Total read pairs with mutations:" in line1:
                    n_depth_1 = int(line1.split(":")[1])
                    n_depth_2 = int(line2.split(":")[1])
                    total_depth = n_depth_1 + n_depth_2
                    merged_header += f"#Total read pairs with mutations:{total_depth}\n"
                    continue
                merged_header += line1.strip() + "; " + line2
        output_f = open(os.path.join(output, f"counts_sample_{row['Sample ID_run1']}-{row['Sample ID_run2']}.csv"), "w")
        output_f.write(merged_header)
        # skip header and read the mutations
        df1 = pd.read_csv(mut_count_f1, skiprows=18)
        df2 = pd.read_csv(mut_count_f2, skiprows=18)
        # for each file find the corresponded wt file and read the wt
        merge = [df1, df2]
        merge = pd.concat(merge)
        merge_counts = merge.groupby(by=["HGVS"], as_index=False)["count"].sum().sort_values("count", ascending=False)
        # save header and counts to file
        merge_counts.to_csv(output_f, mode="a", index=False)
        output_f.close()
        
        cov_f1 = os.path.join(input_dir1, f"coverage_{row['Sample ID_run1']}.csv")
        cov_f2 = os.path.join(input_dir2, f"coverage_{row['Sample ID_run2']}.csv")

        # read coverage files
        # sum up two df
        cov_d1 = pd.read_csv(cov_f1).set_index("pos")
        cov_d2 = pd.read_csv(cov_f2).set_index("pos")
        df_sum = cov_d1.add(cov_d2, fill_value=0)
        output_cov = os.path.join(output, f"coverage_{row['Sample ID_run1']}-{row['Sample ID_run2']}.csv")
        df_sum.to_csv(output_cov)
    return new_samples


def read_json(json_dict):
    """
    Read json input and return data frame of samples and their conditions
    """
    samples = json_dict["samples"]
    condition_def = pd.DataFrame(json_dict["conditions"]["definitions"])
    # get nonselect conditions
    # get wt for nonselect
    try:
        select_relation = condition_def.loc[condition_def["Relationship"].str.contains("is_selection_for")]
        wt_relation = condition_def.loc[condition_def["Relationship"].str.contains("is_wt_control_for")]
        select_relation = select_relation[["Condition 1", "Condition 2"]]
        wt_relation = wt_relation[["Condition 1", "Condition 2"]]
        relation = pd.merge(select_relation, wt_relation, how="left", left_on="Condition 2", right_on="Condition 2")
        relation.columns = ["S", "NS", "WT"]
    except:
        relation = pd.DataFrame({})
    samples_df = pd.DataFrame(samples)
    return relation, samples_df


def parse_params(param1, param2, output, main_log):
    """
    Given two parameter sheets, find samples in common and output
    df for later merging
    :param param1: parameter sheet for the first run
    :param param2: parameter sheet for the second run
    :param main_log: logger
    :return map_df: df with two cols [samples_1, samples_2]
    """
    with open(param1, "r") as p1:
        with open(param2, "r") as p2:
            param_one = json.load(p1)
            param_two = json.load(p2)
            wt_relation1, samples1 = read_json(param_one)
            wt_relation2, samples2 = read_json(param_two)

    # from these two files find samples in the same tile, with the same condition names
    # same replicates and tp
    # log how many samples will be merged
    merge = pd.merge(samples1, samples2, how="inner", on=["Tile ID", "Condition", "Time point", "Replicate"],
                     suffixes=["_run1", "_run2"])
    unique_Tid = merge["Tile ID"].drop_duplicates().tolist()
    unique_Conditions = merge["Condition"].drop_duplicates().tolist()
    merge["new_samples"] = merge["Sample ID_run1"] + "-" + merge["Sample ID_run2"]
    # save new sample names to output
    merge[["new_samples", "Tile ID", "Condition", "Time point", "Replicate"]].to_csv(os.path.join(output, "new_sample_names.csv"), index=False)
    main_log.info(f"{merge.shape[0]} samples will be merged")
    main_log.info(f"Tiles: {unique_Tid}")
    main_log.info(f"Condition: {unique_Conditions}")

    merge_wt = pd.concat([wt_relation1, wt_relation2]).drop_duplicates()
    return merge, merge_wt


def subtract_wt(input_dir1, input_dir2, map_df, wt_map, output):
    """

    """
    for index, row in map_df.iterrows():
        mut_count_f1 = os.path.join(input_dir1, f"counts_sample_{row['Sample ID_run1']}.csv")
        mut_count_f2 = os.path.join(input_dir2, f"counts_sample_{row['Sample ID_run2']}.csv")
        # process select, nonselect and wt together for each sample
        if row["Condition"] in wt_map["S"].tolist():
            print(row["Condition"])
            ns_name = wt_map[wt_map["S"] == row["Condition"]]["NS"].tolist()[0]
            wt_name = wt_map[wt_map["S"] == row["Condition"]]["WT"].tolist()[0]
            # corresponded wt and non-select
            wt = map_df.loc[(map_df["Tile ID"] == row["Tile ID"]) & (map_df["Time point"] == row["Time point"]) & (map_df["Replicate"] == row["Replicate"])]
            wt_run1_sample = wt[wt["Condition"] == wt_name]["Sample ID_run1"].tolist()[0]
            wt_run2_sample = wt[wt["Condition"] == wt_name]["Sample ID_run2"].tolist()[0]
            mut_count_wt1 = os.path.join(input_dir1, f"counts_sample_{wt_run1_sample}.csv")
            mut_count_wt2 = os.path.join(input_dir2, f"counts_sample_{wt_run2_sample}.csv")

            ns = map_df.loc[(map_df["Tile ID"] == row["Tile ID"]) & (map_df["Time point"] == row["Time point"]) & (
                        map_df["Replicate"] == row["Replicate"])]
            ns_run1_sample = ns[ns["Condition"] == ns_name]["Sample ID_run1"].tolist()[0]
            ns_run2_sample = ns[ns["Condition"] == ns_name]["Sample ID_run2"].tolist()[0]
            mut_count_ns1 = os.path.join(input_dir1, f"counts_sample_{ns_run1_sample}.csv")
            mut_count_ns2 = os.path.join(input_dir2, f"counts_sample_{ns_run2_sample}.csv")

            # remove wt counts from run 1
            updated_counts_run1, updated_counts_run1_ns,  updated_wt_run1 = remove_wt_counts(mut_count_f1,
                                                                                            mut_count_wt1, mut_count_ns1)
            updated_counts_run2, updated_counts_run2_ns, updated_wt_run2 = remove_wt_counts(mut_count_f2,
                                                                                            mut_count_wt2, mut_count_ns2)

            # find corresponded merged files and append df to the file
            mut_count_file = os.path.join(output, f"counts_sample_{row['Sample ID_run1']}-{row['Sample ID_run2']}.csv")
            mut_count_wt_file = os.path.join(output, f"counts_sample_{wt_run1_sample}-{wt_run2_sample}.csv")
            mut_count_ns_file = os.path.join(output, f"counts_sample_{ns_run1_sample}-{ns_run2_sample}.csv")

            # write mutations to selection condition
            merge = [updated_counts_run1, updated_counts_run2]
            merge = pd.concat(merge)
            merge_counts = merge.groupby(by=["HGVS"], as_index=False)["count"].sum().sort_values("count",
                                                                                                 ascending=False)
            merge_counts["count"] = merge_counts["count"].astype(int)
            merge_counts.to_csv(mut_count_file, mode="a", index=False)

            # write mutations to wt condition
            merge = [updated_wt_run1, updated_wt_run2]
            merge = pd.concat(merge)
            merge_counts = merge.groupby(by=["HGVS"], as_index=False)["count"].sum().sort_values("count",
                                                                                                 ascending=False)
            merge_counts["count"] = merge_counts["count"].astype(int)
            merge_counts.to_csv(mut_count_wt_file, mode="a", index=False)

            # write mutations to ns condition
            merge = [updated_counts_run1_ns, updated_counts_run2_ns]
            merge = pd.concat(merge)
            merge_counts = merge.groupby(by=["HGVS"], as_index=False)["count"].sum().sort_values("count",
                                                                                                 ascending=False)
            merge_counts["count"] = merge_counts["count"].astype(int)
            merge_counts.to_csv(mut_count_ns_file, mode="a", index=False)


def write_header(input_dir1, input_dir2, samples_df, output):
    """

    """
    new_samples = []
    for index, row in samples_df.iterrows():
        merged_header = ""
        mut_count_f1 = os.path.join(input_dir1, f"counts_sample_{row['Sample ID_run1']}.csv")
        mut_count_f2 = os.path.join(input_dir2, f"counts_sample_{row['Sample ID_run2']}.csv")
        new_samples.append(f"{row['Sample ID_run1']}-{row['Sample ID_run2']}")
        # join headers
        with open(mut_count_f1) as file1, open(mut_count_f2) as file2:
            for line1, line2 in zip(file1, file2):
                if "#" not in line1: break
                if "Final read-depth:" in line1:
                    n_depth_1 = int(line1.split(":")[1])
                    n_depth_2 = int(line2.split(":")[1])
                    total_depth = n_depth_1 + n_depth_2
                    merged_header += f"#Final read-depth:{total_depth}\n"
                    continue
                if "Number of read pairs without mutations:" in line1:
                    n_depth_1 = int(line1.split(":")[1])
                    n_depth_2 = int(line2.split(":")[1])
                    total_depth = n_depth_1 + n_depth_2
                    merged_header += f"#Number of read pairs without mutations:{total_depth}\n"
                    continue
                if "Total read pairs with mutations:" in line1:
                    n_depth_1 = int(line1.split(":")[1])
                    n_depth_2 = int(line2.split(":")[1])
                    total_depth = n_depth_1 + n_depth_2
                    merged_header += f"#Total read pairs with mutations:{total_depth}\n"
                    continue
                merged_header += line1.strip() + "; " + line2
        output_f = open(os.path.join(output, f"counts_sample_{row['Sample ID_run1']}-{row['Sample ID_run2']}.csv"), "w")
        output_f.write(merged_header)
        output_f.close()


def remove_wt_counts(mut_count_f, wt_f, ns_f):
    """

    """
    # load  df
    # skip header and read the mutations
    df = pd.read_csv(mut_count_f, skiprows=18)
    # load wt df
    # skip header and read the mutations
    wt_df = pd.read_csv(wt_f, skiprows=18)
    # load ns df
    # skip header and read the mutations
    ns_df = pd.read_csv(ns_f, skiprows=18)

    # from non-select, remove wt counts
    merged_df_ns = pd.merge(ns_df, wt_df, on="HGVS", how="outer", suffixes=["_ns", "_wt"])
    merged_df_ns = merged_df_ns.fillna(0)
    # subtract wt
    merged_df_ns["updated_ns"] = merged_df_ns["count_ns"] - merged_df_ns["count_wt"]

    # for all the values <= 0, assign 0 to the updated wt
    merged_df_ns["updated_wt"] = merged_df_ns["count_wt"]

    merged_df_ns.loc[(merged_df_ns["updated_ns"] <= 0) & (merged_df_ns["count_ns"] != 0), "updated_wt"] = 0
    merged_df_ns.loc[(merged_df_ns["updated_ns"] <= merged_df_ns["count_ns"]) & (merged_df_ns["updated_ns"] > 0),
                     "updated_wt"] = 0

    merged_df_ns.loc[(merged_df_ns["updated_ns"] <= 0) & (merged_df_ns["count_ns"] > 0), "updated_ns"] = 0
    merged_df_ns.loc[merged_df_ns["count_ns"] == 0, "updated_ns"] = 0
    updated_counts_ns = merged_df_ns.loc[merged_df_ns["updated_ns"] != 0][["HGVS", "updated_ns"]]
    updated_counts_ns.columns = ["HGVS", "count"]
    updated_counts_wt_ns = merged_df_ns.loc[merged_df_ns["updated_wt"] != 0][["HGVS", "updated_wt"]]
    updated_counts_wt_ns.columns = ["HGVS", "ns_count"]

    # from select, remove wt counts
    merged_df = pd.merge(df, wt_df, on="HGVS", how="outer", suffixes=["_s", "_wt"])
    merged_df = merged_df.fillna(0)
    # subtract wt
    merged_df["updated_s"] = merged_df["count_s"] - merged_df["count_wt"]

    # for all the values <= 0, assign 0 to the updated wt
    merged_df["updated_wt"] = merged_df["count_wt"]

    merged_df.loc[(merged_df["updated_s"] <= 0) & (merged_df["count_s"] != 0), "updated_wt"] = 0
    merged_df.loc[(merged_df["updated_s"] <= merged_df["count_s"]) & (merged_df["updated_s"] > 0), "updated_wt"] = 0

    merged_df.loc[(merged_df["updated_s"] <= 0) & (merged_df["count_s"] > 0), "updated_s"] = 0
    merged_df.loc[merged_df["count_s"] == 0, "updated_s"] = 0
    updated_counts = merged_df.loc[merged_df["updated_s"] != 0][["HGVS", "updated_s"]]
    updated_counts.columns = ["HGVS", "count"]

    updated_counts_wt_s = merged_df.loc[merged_df["updated_wt"] != 0][["HGVS", "updated_wt"]]
    updated_counts_wt_s.columns = ["HGVS", "count"]

    # join updated counts wt ns and s into new wt
    # the updated wt counts contains the variants that are not found in ns or s
    updated_counts_wt = pd.merge(updated_counts_wt_ns, updated_counts_wt_s, how="inner", on="HGVS")
    return updated_counts, updated_counts_ns, updated_counts_wt


def merge_main(args):

    # check args and make corresponding output dir
    # log input/output information
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    logging = logginginit(args.log_level, os.path.join(args.output, "mergeRuns.log"))
    main_log = logging.getLogger("mergeRuns.log")
    main_log.info(f"Input dir 1: {args.dir1} with Parameter Sheet: {args.paramOne}")
    main_log.info(f"Input dir 2: {args.dir2} with Parameter Sheet: {args.paramTwo}")
    main_log.info(f"Merges mutation count file will be saved to {args.output}")
    map_df, wt_map_df = parse_params(args.paramOne, args.paramTwo, args.output, main_log)
    if args.subtractWT:
        main_log.info("subtract WT turned on, wt counts will be subtracted from select/non-select conditions")
        write_header(args.dir1, args.dir2, map_df, args.output)
        subtract_wt(args.dir1, args.dir2, map_df, wt_map_df, args.output)

    else:
        # map_df = parse_params(args.paramOne, args.paramTwo, args.output, main_log)
        new_samples = merge_runs(args.dir1, args.dir2, map_df, args.output)

    main_log.info(f"Merged files are saved to {args.output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TileSeq mutation counts')
    # user input arguments
    parser.add_argument("-p1", "--paramOne", help="Path to parameter sheet (JSON) for the first run", type=str, required=True)
    parser.add_argument("-p2", "--paramTwo", help="Path to parameter sheet (JSON) for the second run", type=str, required=True)
    parser.add_argument("-d1", "--dir1", help="Path to *_mut_count folder for the first run", type=str, required=True)
    parser.add_argument("-d2", "--dir2", help="Path to *_mut_count folder for the second run", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output folder path (this folder will be made if not exists when "
                                               "running the script",type=str, required=True)
    parser.add_argument("--subtractWT", help="Subtract wt counts from non-select before merging", action="store_true")
    parser.add_argument("-log", "--log_level", help="set log level: debug, \
            info, warning, error, critical. (default = debug)", type=str, default="debug")
    # test dirs (local, CHEK2)
    # dir r2r3 /Users/roujia/Desktop/Tileseq_DMS/CHEK2/210216_CHEK2_R2R3/210216_CHEK2_R2R3_2021-02-23-18-39-16_mut_count
    # dir 200310 /Users/roujia/Desktop/Tileseq_DMS/CHEK2/200310_rerun_CHEK2/v419_2021-02-02-17-57-49_mut_count
    # param1 210216_CHEK2.json
    # param2 200310_CHECK2_newtile11.json

    args = parser.parse_args()

    merge_main(args)
