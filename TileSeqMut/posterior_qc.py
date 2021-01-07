#!/usr/bin/env python3.7

import sys
import os
import pandas as pd
import argparse
from ast import literal_eval
from matplotlib import rc

import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
sys.path.append('..')
# from TileSeqMut import help_functions
import help_functions

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

# QC script for posteriors


class PosteriorQC(object):

    def __init__(self, mut_count_dir, param_json, log_level="INFO"):
        """
        @param mut_count_dir: Directory ends with _mut_count
        @param param_json: Json parameter sheet
        """
        # validate input parameters
        if not os.path.isdir(mut_count_dir):
            raise NotADirectoryError(f"Dir does not exist: {mut_count_dir}")

        if not os.path.isfile(param_json):
            raise FileNotFoundError(f"File does not exist: {param_json}")

        self._project, self._seq, self._cds_seq, self._tile_map, \
            self._region_map, self._samples, self._var = help_functions.parse_json(param_json)

        self._mut_dir = mut_count_dir
        self._output = os.path.join(mut_count_dir, "posterior_qc")

        if not os.path.isdir(self._output):
            os.mkdir(self._output)

        log_object = help_functions.logginginit(log_level, "posteriorQC.log")
        self._logger = log_object.getLogger("posteriorQC")

        self._logger.info(f"Reading files from {os.path.abspath(self._mut_dir)}")
        self._logger.info(f"Reading parameters from {os.path.abspath(param_json)}")
        self._logger.info(f"Saving output plots to {os.path.abspath(self._output)}")

    def prob_analysis(self, row):
        """
        @param row: row of dataframe contains path to prob files and sample info
        """
        # load prob_df
        fig_title = f"{row['Condition']}_tile{row['Tile ID']}_rep{row['Replicate']}_t{row['Time point']}"
        prob_df = pd.read_csv(row.prob, names=["mut", "prob", "read"], header=None)
        # print(prob_df)
        # self._plot_prob(prob_df)
        # self._plot_marginal(prob_df)
        all_prob_df = pd.read_csv(row.prob_all, names=["mut", "prob", "read"], header=None)
        self._plot_prob(all_prob_df, fig_title)
        # plot_name = f'Tile{row["Tile ID"]}_{row["Condition"]}_t{row["Time point"]}_rep{row["Replicate"]}.png'
        # self._plot_counts(prob_df, plot_name)
        # plot_name = f'Tile{row["Tile ID"]}_{row["Condition"]}_t{row["Time point"]}_rep{row["Replicate"]}_all.png'
        # self._plot_counts(prob_df, plot_name)

    def _plot_prob(self, df, fig_title):
        """
        Plot probability of read 1 and read 2
        """
        print(fig_title)
        # mutations only on read 1
        read1_only = df.loc[df["read"] == "r1"]
        read1_only_marginal = self._calculate_frequency(read1_only)
        # mutations only on read 2
        read2_only = df.loc[df["read"] == "r2"]
        read2_only_marginal = self._calculate_frequency(read2_only)
        print(read2_only_marginal)
        print(read1_only_marginal)
        # plot marginal
        sns.scatterplot(read1_only_marginal.index, read1_only_marginal.marginal_freq, s = 8, label="Read 1")
        sns.scatterplot(read2_only_marginal.index, read2_only_marginal.marginal_freq, s = 8, label="Read 2")
        plt.title(fig_title)
        plt.savefig(os.path.join(self._output, f"{fig_title}_marginal_r1r2.png"))
        plt.close()
        # mutations found on both reads
        both_read = df.loc[(df["read"] != "r1") & (df["read"] != "r2")]
        both_read['read'] = [literal_eval(x) for x in both_read['read']]
        both_read['prob'] = [literal_eval(x) for x in both_read['prob']]
        both_read[["prob_r1", "prob_r2"]] = pd.DataFrame(both_read['prob'].tolist(), index=both_read.index)
        both_read[["label_r1", "label_r2"]] = pd.DataFrame(both_read['read'].tolist(), index=both_read.index)

        read1_only.prob = read1_only.prob.astype(float)
        read2_only.prob = read2_only.prob.astype(float)
        read1_only = read1_only.sort_values(by="prob", ascending=False).reset_index()
        read2_only = read2_only.sort_values(by="prob", ascending=False).reset_index()

        read1_only["perc"] = read1_only.index / (read1_only.shape[0] - 1)
        read2_only["perc"] = read2_only.index / (read2_only.shape[0] - 1)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,7))
        sns.scatterplot(read1_only.perc, read1_only.prob, label="Read 1", ax=ax1, s=7)
        sns.scatterplot(read2_only.perc, read2_only.prob, label="Read 2", ax=ax1, s=7)

        # plot prob correlation of r1 and r2
        sns.scatterplot(both_read.prob_r1, both_read.prob_r2, s=7, ax=ax2, color="black")
        pcc = stats.pearsonr(both_read.prob_r1, both_read.prob_r2)[0]
        ax1.set_title("Posterior prob of mutations found only on Read 1/Read 2")
        ax1.set_ylabel("Posterior prob", fontsize=12)
        ax1.set_xlabel("Perc. of mutations", fontsize=12)

        ax2.set_title("Corr of posterior prob of mutations found on both reads")
        ax2.text(0, 1, f"PCC: {float('{:.4f}'.format(pcc))}")
        ax2.set_ylabel("Posterior prob (Read 1)", fontsize=12)
        ax2.set_xlabel("Posterior prob (Read 2)", fontsize=12)

        fig.suptitle(fig_title, fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(self._output, f"{fig_title}.png"))
        plt.close()

        # # make plot for mutation on r1 only and r2 only vs all
        # n_r1 = read1_only.mut.value_counts().to_frame().reset_index()
        # n_r2 = read2_only.mut.value_counts().to_frame().reset_index()
        # print(n_r1)
        # print(n_r2)
        #
        # sns.scatterplot(n_r1.index, n_r1["mut"], label="Read 1")
        # sns.scatterplot(n_r2.index, n_r2["mut"], label="Read 2")
        # # bar_df = pd.DataFrame({"type": ["Mutations on read1", "Mutations on read2", "Mutations on both reads"],
        # #                        "count": [n_r1, n_r2, both]})
        # # sns.barplot(x="type", y="count", data=bar_df, palette=["#F07D78", "#B7F86E", "#316FA3"])
        # # plt.title("Number of mutations (before filtering) on each read")
        # plt.show()
        # plt.close()
        # exit()

    def _plot_marginal(self, df):
        """

        """
        # mutations only on read 1
        read1_only = df.loc[df["read"] == "r1"]
        # mutations only on read 2
        read2_only = df.loc[df["read"] == "r2"]

        # mutations found on both reads
        both_read = df.loc[(df["read"] != "r1") & (df["read"] != "r2")]
        both_read['read'] = [literal_eval(x) for x in both_read['read']]
        both_read['prob'] = [literal_eval(x) for x in both_read['prob']]
        both_read[["prob_r1", "prob_r2"]] = pd.DataFrame(both_read['prob'].tolist(), index=both_read.index)
        both_read[["label_r1", "label_r2"]] = pd.DataFrame(both_read['read'].tolist(), index=both_read.index)

        r1_marginal = self._calculate_frequency(read1_only)
        r1_marginal["read"] = "r1"
        r2_marginal = self._calculate_frequency(read2_only)
        r2_marginal["read"] = "r2"
        both_read_marginal = self._calculate_frequency(both_read)
        both_read_marginal["read"] = "both"

        marginal_concat = pd.concat([r1_marginal, r2_marginal, both_read_marginal])
        sns.violinplot(x = "read", y="marginal_freq", data=marginal_concat)
        plt.show()
        plt.close()
        exit()

    def _calculate_frequency(self, df):
        """
        Calculate marginal frequencies for all the variants in the table
        """
        # count each unique variant in the table
        mut_counts = df["mut"].value_counts()
        mut_counts = mut_counts.to_frame().reset_index()
        mut_counts.columns = ["mut", "count"]
        total_mut = mut_counts["count"].sum()
        mut_counts["marginal_freq"] = mut_counts["count"] / total_mut
        return mut_counts

    def _plot_counts(self, df, f_name):
        """
        make stacked barplots of number of mutations on each read
        """
        # mutations only on read 1
        read1_only = df.loc[df["read"] == "r1"]
        # mutations only on read 2
        read2_only = df.loc[df["read"] == "r2"]

        # mutations found on both reads
        both_read = df.loc[(df["read"] != "r1") & (df["read"] != "r2")]
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
        plt.savefig(os.path.join(self._output, f_name))
        plt.close()

    def runQC(self):
        """
        Run posterior QC, generate posterior QC plots in a subdir
        """
        # go through sample ID list and find corresponding files
        f_table = []
        # join samples table into one col with _ delim
        for sample in self._samples["Sample ID"].tolist():
            # filtered prob df
            filtered_prob_f = os.path.join(self._mut_dir, f"{sample}_posprob.csv")
            if not os.path.isfile(filtered_prob_f):
                raise FileNotFoundError(f"{sample}_posprob.csv")
            all_prob_f = os.path.join(self._mut_dir, f"{sample}_posprob_all.csv")
            if not os.path.isfile(all_prob_f):
                raise FileNotFoundError(f"{sample}_posprob_all.csv")
            f_table.append([sample, filtered_prob_f, all_prob_f])

        f_df = pd.DataFrame(f_table, columns=["Sample ID", "prob", "prob_all"])
        # merge this to sample df
        merge_file_samples = pd.merge(self._samples, f_df, on="Sample ID", how="left")

        # go through each condition to make QC plots
        conditions = list(set(merge_file_samples["Condition"]))
        for c in conditions:
            print(c)
            # get number of subplots we have to make for this condition
            # fig, axes = plt.subplots(4, 2, figsize=(13, 30))
            c_df = merge_file_samples.loc[merge_file_samples["Condition"] == c]
            # sort df by tiles so the output plots are in the right order
            c_df = c_df.sort_values(by="Tile ID")
            # print(c_df["Sample ID"])
            # exit()
            c_df.apply(lambda x: self.prob_analysis(x), axis=1)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TileSeq mutation posterior QC')
    # user input arguments
    parser.add_argument("-i", "--input", help="Path to folder that contains mutation counts (*_mut_count)", type=str,
                        required=True)
    parser.add_argument("-p", "--param", help="json paramter file", type=str, required=True)

    args = parser.parse_args()

    p_qc = PosteriorQC(args.input, args.param)
    p_qc.runQC()
