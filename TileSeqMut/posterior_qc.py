#!/usr/bin/env python3.7
import glob
import math
import sys
import os
import pandas as pd
import numpy as np
import argparse
from ast import literal_eval
from matplotlib_venn import venn2
from fpdf import FPDF
import progressbar

from matplotlib import rc

import seaborn as sns
import warnings
from scipy import stats
import matplotlib.pyplot as plt
sys.path.append('..')
from TileSeqMut import help_functions
# import help_functions
warnings.simplefilter(action='ignore')
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
        self._i = 0

    def prob_analysis(self, row, pdf, bar):
        """
        @param row: row of dataframe contains path to prob files and sample info
        """
        # load prob_df
        fig_title = f"{row['Condition']}_tile{row['Tile ID']}_rep{row['Replicate']}_t{row['Time point']}"
        # prob_df = pd.read_csv(row.prob, names=["mut", "prob", "read"], header=None)
        # print(prob_df)
        # self._plot_prob(prob_df)
        # self._plot_marginal(prob_df)
        sns.set(font_scale=1.7)
        sns.set_style("whitegrid")
        fig, axes = plt.subplots(3, 2, figsize=(16, 18))
        axes[-1, -1].axis('off')
        all_prob_df = pd.read_csv(row.prob_all, names=["mut", "prob", "read", "pass"], header=None)
        pop_r1_df = pd.read_csv(row.popmut_r1, skiprows=9)
        pop_r2_df = pd.read_csv(row.popmut_r2, skiprows=9)
        self._plot_prob(all_prob_df, fig_title, fig, axes)
        self._plot_r1r2_popmut(pop_r1_df, pop_r2_df, fig_title, fig, axes)
        # plt.tight_layout(pad=0.6)
        fig.suptitle(fig_title, fontweight='bold')
        fig.tight_layout(pad=1)
        plt.savefig(os.path.join(self._output, f"{fig_title}.png"), format="png")
        plt.close()
        self._logger.debug("plot closed")
        pdf.add_page()
        pdf.image(os.path.join(self._output, f"{fig_title}.png"), 2, 25, 200, 235)
        bar.update(self._i + 1)
        self._logger.debug("added to pdf")
        self._i +=1
        # self._logger.info(f"Processed sample {fig_title}")

        # plot_name = f'Tile{row["Tile ID"]}_{row["Condition"]}_t{row["Time point"]}_rep{row["Replicate"]}.png'
        # self._plot_counts(prob_df, plot_name)
        # plot_name = f'Tile{row["Tile ID"]}_{row["Condition"]}_t{row["Time point"]}_rep{row["Replicate"]}_all.png'
        # self._plot_counts(prob_df, plot_name)

    def _plot_prob(self, df, fig_title, fig, axes):
        """
        Plot probability of read 1 and read 2
        """

        # mutations only on read 1
        read1_only = df.loc[df["read"] == "r1"]

        # mutations only on read 2
        read2_only = df.loc[df["read"] == "r2"]

        read1_only.prob = read1_only.prob.astype(float)
        read2_only.prob = read2_only.prob.astype(float)
        # self._logger.info("Calculating frequencies..")
        read1_only_marginal = self._calculate_frequency(read1_only)
        read2_only_marginal = self._calculate_frequency(read2_only)
        # self._logger.info("Finished calc frequencies..")
        map_r1 = {"1": "R1_Passed", "-1": "R1_Discarded"}
        map_r2 = {"1": "R2_Passed", "-1": "R2_Discarded"}

        read1_only = read1_only.replace({"pass": map_r1})
        read2_only = read2_only.replace({"pass": map_r2})
        # self._logger.info("Finished mapping values..")
        top90_cut = np.log(0.9)
        top50_cut = np.log(0.5)
        # plot counts vs avg prob

        sns.scatterplot(read1_only_marginal["marginal_freq"], np.log(read1_only_marginal["prob"]),
                        label="Read 1", ax=axes[0][1], alpha=0.7)
        sns.scatterplot(read2_only_marginal["marginal_freq"], np.log(read2_only_marginal["prob"]),
                        label="Read 2", ax=axes[0][1], alpha=0.7)
        max_x = max(read1_only_marginal["marginal_freq"].max(), read2_only_marginal["marginal_freq"].max())
        # max_x = math.floor(max_x*10)/10

        axes[0][1].axhline(top90_cut, ls='--', alpha=.7, color="#FA5D5D")
        axes[0][1].axhline(top50_cut, ls='--', alpha=.7, color="#FA5D5D")

        axes[0][1].text((max_x/3)*2, top50_cut-0.38, "posterior = 0.5", fontsize=10, color="grey")
        axes[0][1].text((max_x/3)*2, top90_cut-0.38, "posterior = 0.9", fontsize=10, color="grey")

        axes[0][1].set_ylabel("Average posterior prob. (log)")
        axes[0][1].set_xlabel("Marginal frequency")
        axes[0][1].set_title("mutations found only on R1 or R2")

        # axes[0][0].get_legend().remove()
        # plt.show()

        # plt.close()
        # plot marginal
        # sns.scatterplot(read1_only_marginal.index, read1_only_marginal.marginal_freq, s=20, ax=axes[0][0],
        #                 label="Read 1", linewidth=0.000001, alpha = 0.6)
        # sns.scatterplot(read2_only_marginal.index, read2_only_marginal.marginal_freq, s=20, ax=axes[0][0],
        #                 label="Read 2", linewidth=0.000001, alpha = 0.6)
        #
        # axes[0][0].set_title("Marginal frequencies of mutations found only on R1 or R2")
        # axes[0][0].set_xlabel("Unique variants")
        # axes[0][0].set_ylabel("Marginal Frequency")

        # plt.title(fig_title)
        # plt.savefig(os.path.join(self._output, f"{fig_title}_marginal_r1r2.pdf"))
        # plt.close()
        # mutations found on both reads
        # self._logger.info("Reading both read mutations..")
        both_read = df.loc[(df["read"] != "r1") & (df["read"] != "r2")]
        both_read['read'] = [literal_eval(x) for x in both_read['read']]
        both_read['prob'] = [literal_eval(x) for x in both_read['prob']]

        both_read[["prob_r1", "prob_r2"]] = pd.DataFrame(both_read['prob'].tolist(), index=both_read.index)
        both_read[["label_r1", "label_r2"]] = pd.DataFrame(both_read['read'].tolist(), index=both_read.index)
        both_read["pass"] = "Discarded"
        both_read["label"] = "both"
        both_read.loc[(both_read.prob_r1 > 0.5) & (both_read.prob_r2 > 0.5), "pass"] = "Passed"
        read1_only["log_prob"] = np.log(read1_only.prob)
        read2_only["log_prob"] = np.log(read2_only.prob)
        # self._logger.info("started sorting...")
        read1_only = read1_only.sort_values(by="prob", ascending=False).reset_index()
        read2_only = read2_only.sort_values(by="prob", ascending=False).reset_index()
        # self._logger.info("finished sorting...")
        read1_only["label"] = "r1_only"
        read2_only["label"] = "r2_only"

        read1_only["perc"] = read1_only.index / (read1_only.shape[0] - 1)
        read2_only["perc"] = read2_only.index / (read2_only.shape[0] - 1)
        # unique variants r1 only
        unique_r1_only = read1_only.drop_duplicates(subset=["mut"])
        # unique variants r2 only
        unique_r2_only = read2_only.drop_duplicates(subset=["mut"])
        # unique variants both
        unique_both = both_read.drop_duplicates(subset=["mut"])

        unique_r1_only[['pos', 'ref', 'alt', 'qual']] = unique_r1_only['mut'].str.split('|', expand=True)
        unique_r2_only[['pos', 'ref', 'alt', 'qual']] = unique_r2_only['mut'].str.split('|', expand=True)
        unique_both[['pos', 'ref', 'alt', 'qual']] = unique_both['mut'].str.split('|', expand=True)

        # print(unique_r1_only.columns)
        # unique variants r1 only
        unique_r1_only = unique_r1_only.drop_duplicates(subset=['pos', 'ref', 'alt'])
        # unique variants r2 only
        unique_r2_only = unique_r2_only.drop_duplicates(subset=['pos', 'ref', 'alt'])
        # unique variants both
        unique_both = unique_both.drop_duplicates(subset=['pos', 'ref', 'alt'])

        unique_r1_only["type"] = "SNP"
        unique_r1_only.loc[(unique_r1_only["alt"].str.contains("ins")) | (unique_r1_only["alt"].str.contains("ins")),
                        "type"] = "INDEL"

        unique_r2_only["type"] = "SNP"
        unique_r2_only.loc[(unique_r2_only["alt"].str.contains("ins")) | (unique_r2_only["alt"].str.contains(
            "ins")), "type"] = "INDEL"

        unique_both["type"] = "SNP"
        unique_both.loc[(unique_r1_only["alt"].str.contains("ins")) | (unique_both["alt"].str.contains(
            "ins")), "type"] = "INDEL"

        # self._logger.info("value counts...")
        total_both = unique_both[["label", "pass"]].value_counts().to_frame()
        total_r1 = unique_r1_only[["label", "pass"]].value_counts().to_frame()
        total_r2 = unique_r2_only[["label", "pass"]].value_counts().to_frame()
        bar_plot_df = pd.concat([total_both, total_r2, total_r1]).reset_index()
        bar_plot_df.columns = ["label", "pass", "count"]
        # self._logger.info("Finished value counts...")

        map = {"R1_Passed": "Passed", "R1_Discarded": "Discarded", "R2_Passed": "Passed", "R2_Discarded": "Discarded"}
        bar_plot_df = bar_plot_df.replace({"pass": map})
        color_dict = {"Passed": "#45736A", "Discarded": "#D94E4E"}

        g = sns.barplot(x="label", y="count", hue="pass", data=bar_plot_df, ax=axes[0][0], palette=color_dict)
        for p in g.patches:
            height = p.get_height()
            if math.isnan(height):
                height=0
            g.text(p.get_x() + p.get_width() / 2.,
                    height + 3,
                    int(height),
                    ha="center")
        handles, labels = axes[0][0].get_legend_handles_labels()
        axes[0][0].legend(handles=handles, labels=labels)
        axes[0][0].set_xlabel("")
        axes[0][0].set_title("Number of unique mutations found")

        # self._logger.info("Join R1 and R2 2/3nt...")
        # join r1 and r2
        join_R1R2 = [read1_only, read2_only]
        join_R1R2_df = pd.concat(join_R1R2)
        # self._logger.info("Joined R1 and R2 2/3nt...")
        # Create an array with the colors you want to use
        colors = ["#6BA3F2", "#3866A6"]
        # Set your custom color palette
        customPalette = sns.set_palette(sns.color_palette(colors))
        #
        colors = ["#F2C46B", "#A67B28"]
        customPalette2 = sns.set_palette(sns.color_palette(colors))
        color_dict = {"R1_Passed": "#6BA3F2", "R1_Discarded": "#3866A6", "R2_Passed": "#F2C46B",
                      "R2_Discarded":"#A67B28"}
        sns.scatterplot(x="perc", y="log_prob", hue="pass", palette=color_dict, data=join_R1R2_df,
                        ax=axes[1][0], s=20, style="pass", linewidth=0.000001, alpha = 0.7)
        handles, labels = axes[1][0].get_legend_handles_labels()
        axes[1][0].legend(handles=handles, labels=labels)
        # sns.scatterplot(x="perc", y="log_prob", hue="pass", palette=color_dict, data=read2_only,
        #                 ax=axes[1][0], s=20, style="pass", linewidth=0.000001, alpha = 0.7)
        # sns.scatterplot(x="perc", y="prob", hue="pass", data=read2_only , label=f"Read 2: n = {read2_only.shape[0]}",
        #                 ax=ax1, s=20)
        axes[1][0].axhline(top90_cut, ls='--', alpha=.7, color="#FA5D5D")
        axes[1][0].axhline(top50_cut, ls='--', alpha=.7, color="#FA5D5D")

        axes[1][0].text(0.6666, top50_cut-0.38, "posterior = 0.5", fontsize=10, color="grey")
        axes[1][0].text(0.6666, top90_cut-0.38, "posterior = 0.9", fontsize=10, color="grey")
        axes[1][0].set_title("Posterior prob of mutations only found on R1 or R2")
        axes[1][0].set_ylabel("Posterior prob (log)")
        axes[1][0].set_xlabel("Cumulative proportion of mutations")

        # self._logger.info("plotted plot 3...")


        # # plot prob correlation of r1 and r2
        sns.scatterplot(np.log(both_read.prob_r1), np.log(both_read.prob_r2), s=10, ax=axes[1][1], color="black")
        pcc = stats.pearsonr(both_read.prob_r1, both_read.prob_r2)[0]
        axes[1][1].set_title("Corr. of posterior prob of mutations found on both reads")
        axes[1][1].text(-3, -1, f"PCC: {float('{:.4f}'.format(pcc))}")
        axes[1][1].set_ylabel("Posterior prob (log) (Read 1)")
        axes[1][1].set_xlabel("Posterior prob (log) (Read 2)")

        # self._logger.info("plotted plot 4...")

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

    def _plot_r1r2_popmut(self, df_r1, df_r2, fig_title, fig, axes):
        """
        For R1 and R2, plot the distribution of popcode mutations based on the margional frequency
        """
        # print(df_r1)
        # print(df_r2)
        if df_r1.empty:
            df_r1 = pd.DataFrame({}, columns=["HGVS", "count"])
        if df_r2.empty:
            df_r2 = pd.DataFrame({}, columns=["HGVS", "count"])
        total_r1 = df_r1["count"].sum()
        df_r1["marginal freq"] = df_r1["count"] / total_r1
        total_r2 = df_r2["count"].sum()
        df_r2["marginal freq"] = df_r2["count"] / total_r2

        # get union set of hgvs
        all_hgvs = set(df_r1["HGVS"].tolist()+df_r2["HGVS"].tolist())
        intersect  = [value for value in df_r1["HGVS"].tolist() if value in df_r2["HGVS"].tolist()]

        unique_r1 = df_r1.loc[~df_r1.HGVS.isin(df_r2)].dropna()
        unique_r2 = df_r2.loc[~df_r2.HGVS.isin(df_r1)].dropna()

        # sort by freq
        unique_r1 = unique_r1.sort_values(by="marginal freq", ascending=False)
        unique_r2 = unique_r2.sort_values(by="marginal freq", ascending=False)

        if len(all_hgvs) != 0:
            j_dist = len(intersect) / len(all_hgvs)
            venn2([set(unique_r1.HGVS.tolist()), set(unique_r2.HGVS.tolist())], set_labels=('Read 1', 'Read 2'),
                  ax=axes[2][0])
            axes[2][0].set_title("Unique 2/3nt changes in R1 oCBSr R2 (passed threshold cutoff)")
            axes[2][0].text(.5, .6, f"Jaccard Index: {round(j_dist, 2)}")

        else: # no 2/3nt changes found on either reads
            axes[2][0].axis("off")


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
        avg_prob = df.groupby(['mut'])['prob'].mean().to_frame().reset_index()
        mut_counts = pd.merge(mut_counts, avg_prob, how="left", on="mut")
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
            # filtered_prob_f = os.path.join(self._mut_dir, f"{sample}_posprob.csv")
            # if not os.path.isfile(filtered_prob_f):
            #     raise FileNotFoundError(f"{sample}_posprob.csv")
            all_prob_f = os.path.join(self._mut_dir, f"{sample}_posprob_all.csv")
            if not os.path.isfile(all_prob_f):
                # raise FileNotFoundError(f"{sample}_posprob_all.csv")
                self._logger(f"Posterior probability files not found for {sample}_posprob_all.csv")
                continue

            pop_r1_f = os.path.join(self._mut_dir, f"counts_sample_{sample}_r1.csv")
            pop_r2_f = os.path.join(self._mut_dir, f"counts_sample_{sample}_r2.csv")
            f_table.append([sample, all_prob_f, pop_r1_f, pop_r2_f])

        f_df = pd.DataFrame(f_table, columns=["Sample ID", "prob_all", "popmut_r1", "popmut_r2"])
        # merge this to sample df
        merge_file_samples = pd.merge(self._samples, f_df, on="Sample ID", how="left")

        # go through each condition to make QC plots
        conditions = list(set(merge_file_samples["Condition"]))
        total_sample = merge_file_samples.shape[0]
        self._logger.info(f"Total sampels = {total_sample}")
        bar = progressbar.ProgressBar(maxval=total_sample, widgets=[progressbar.Bar('=', '[', ']'), ' ',
                                                                    progressbar.Percentage()])
        bar.start()
        for c in conditions:
            # get number of subplots we have to make for this condition
            # fig, axes = plt.subplots(4, 2, figsize=(13, 30))
            c_df = merge_file_samples.loc[merge_file_samples["Condition"] == c]
            # sort df by tiles so the output plots are in the right order
            c_df = c_df.sort_values(by="Tile ID")
            # print(c_df["Sample ID"])
            # exit()
            pdf = FPDF('P', 'mm', 'A4')
            c_df.apply(lambda x: self.prob_analysis(x, pdf, bar), axis=1)
            #
            # # merge all the files for this condition
            # c_list = glob.glob(f"{self._output}/{c}*")
            # c_list.sort()
            # for image in c_list:
            #     pdf.add_page()
            #     pdf.image(image, 2, 25, 200, 245)

            pdf.output(f"{self._output}/{c}_posteriorQC.pdf", "F")

        bar.finish()
        os.system(f"rm {self._output}/*.png")
        self._logger.info("Posterior QC completed.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TileSeq mutation posterior QC')
    # user input arguments
    parser.add_argument("-i", "--input", help="Path to folder that contains mutation counts (*_mut_count)", type=str,
                        required=True)
    parser.add_argument("-p", "--param", help="json paramter file", type=str, required=True)

    args = parser.parse_args()

    p_qc = PosteriorQC(args.input, args.param)
    p_qc.runQC()
