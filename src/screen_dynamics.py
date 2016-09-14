#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import re


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def gRNA_scatter(s1, s2, prefix=""):
    # Scatter of gRNA change
    fig, axis = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8, 8))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2).fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        if "original" in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                b = x["plasmid_pool_TCR"]
            if "WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                b = x["plasmid_pool_WNT"]
        elif "plasmid" not in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                b = x["gDNA_Jurkat"]
            if "_4_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                b = x["gDNA_HEKclone4"]
            if "_6_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                b = x["gDNA_HEKclone6"]
        else:
            b = x[0]

        colors = pd.DataFrame()
        colors[sns.color_palette("colorblind")[0]] = x.index.str.contains("Wnt")
        colors[sns.color_palette("colorblind")[1]] = x.index.str.contains("CTRL")
        colors[sns.color_palette("colorblind")[2]] = x.index.str.contains("Tcr")
        colors[sns.color_palette("colorblind")[3]] = x.index.str.contains("Ess")
        colors = colors.apply(lambda x: x[x].index.tolist()[0], axis=1).tolist()

        lims = [0, np.max([np.log2(1 + x[screen]), np.log2(1 + b)])]
        axis[i].plot(lims, lims, '--', color='black', alpha=0.75)
        axis[i].set_title(screen)

        axis[i].scatter(np.log2(1 + x[screen]), np.log2(1 + b), color=colors, alpha=0.5)

    for i in range(0, len(axis), 2):
        axis[i].set_ylabel("gRNA frequency in plasmid (log2)")
    for ax in axis[-2:]:
        ax.set_xlabel("gRNA frequency in CROP-seq screen (log2)")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.scatter.svg".format(prefix)), bbox_inches="tight")


def gRNA_maplot(s1, s2, prefix=""):
    # Rank of gRNA change
    fig, axis = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8, 8))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2).fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        if "original" in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                M = np.log2(x[screen] * x["plasmid_pool_TCR"]) / 2.
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["plasmid_pool_TCR"])
            if "WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                M = np.log2(x[screen] * x["plasmid_pool_WNT"]) / 2.
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["plasmid_pool_WNT"])
        elif "plasmid" not in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                M = np.log2(x[screen] * x["gDNA_Jurkat"]) / 2.
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_Jurkat"])
            if "_4_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                M = np.log2(x[screen] * x["gDNA_HEKclone4"]) / 2.
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_HEKclone4"])
            if "_6_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                M = np.log2(x[screen] * x["gDNA_HEKclone6"]) / 2.
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_HEKclone6"])
        else:
            M = np.log2((1 + x[screen]) * (1 + x[0])) / 2.
            fc = np.log2(1 + x[screen]) - np.log2(1 + x[0])

        fc.name = screen
        if i == 0:
            xx = pd.DataFrame(fc)
        else:
            xx = xx.join(fc, how="outer")

        colors = pd.DataFrame()
        colors[sns.color_palette("colorblind")[0]] = x.index.str.contains("Wnt")
        colors[sns.color_palette("colorblind")[1]] = x.index.str.contains("CTRL")
        colors[sns.color_palette("colorblind")[2]] = x.index.str.contains("Tcr")
        colors[sns.color_palette("colorblind")[3]] = x.index.str.contains("Ess")
        colors = colors.apply(lambda x: x[x].index.tolist()[0], axis=1).tolist()

        axis[i].scatter(M, fc, color=colors, alpha=0.5)
        axis[i].axhline(y=0, color='black', linestyle='--', lw=0.5)

        axis[i].set_title(screen)

    for i in range(0, len(axis), 2):
        axis[i].set_ylabel("M")
    for ax in axis[-2:]:
        ax.set_xlabel("A")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.maplot.svg".format(prefix)), bbox_inches="tight")


def gRNA_rank(s1, s2, prefix=""):
    # Rank of gRNA change
    fig, axis = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8, 8))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2).fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        if "original" in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["plasmid_pool_TCR"])
            if "WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["plasmid_pool_WNT"])
        elif "plasmid" not in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_Jurkat"])
            if "_4_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_HEKclone4"])
            if "_6_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_HEKclone6"])
        else:
            fc = np.log2(1 + x[screen]) - np.log2(1 + x[0])

        fc.name = screen
        if i == 0:
            xx = pd.DataFrame(fc)
        else:
            xx = xx.join(fc, how="outer")

        colors = pd.DataFrame()
        colors[sns.color_palette("colorblind")[0]] = x.index.str.contains("Wnt")
        colors[sns.color_palette("colorblind")[1]] = x.index.str.contains("CTRL")
        colors[sns.color_palette("colorblind")[2]] = x.index.str.contains("Tcr")
        colors[sns.color_palette("colorblind")[3]] = x.index.str.contains("Ess")
        colors = colors.apply(lambda x: x[x].index.tolist()[0], axis=1).tolist()

        axis[i].scatter(fc.rank(ascending=False, method="first"), fc, color=colors, alpha=0.5)
        axis[i].axhline(y=0, color='black', linestyle='--', lw=0.5)

        axis[i].set_title(screen)

    for i in range(0, len(axis), 2):
        axis[i].set_ylabel("gRNA fold-change")
    for ax in axis[-2:]:
        ax.set_xlabel("gRNA rank")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.svg".format(prefix)), bbox_inches="tight")

    xx.to_csv(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.csv".format(prefix)), index=True)


def gRNA_rank_stimulus(xx, s2, prefix=""):
    # Difference between unstimulated/stimulated
    fig, axis = plt.subplots(1, 3, sharex=False, sharey=True, figsize=(12, 3))
    axis = iter(axis.flatten())

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2).fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        if "original" in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["plasmid_pool_TCR"])
            if "WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["plasmid_pool_WNT"])
        elif "plasmid" not in prefix:
            if "TCR" in screen:
                x = x.ix[x.index[~x.index.str.contains("Wnt")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_Jurkat"])
            if "_4_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_HEKclone4"])
            if "_6_WNT" in screen:
                x = x.ix[x.index[~x.index.str.contains("Tcr")]]
                fc = np.log2(1 + x[screen]) - np.log2(1 + x["gDNA_HEKclone6"])
        else:
            fc = np.log2(1 + x[screen]) - np.log2(1 + x[0])

        fc.name = screen
        if i == 0:
            xx = pd.DataFrame(fc)
        else:
            xx = xx.join(fc, how="outer")

    screens = s2.columns[::-1]
    for i in range(0, len(s2.columns), 2):
        ax = axis.next()
        fc = (xx[screens[i + 1]] - xx[screens[i]]).dropna()

        fc.name = screens[i + 1]
        if i == 0:
            ax.set_ylabel("gRNA fold-change (stimulated / unstimulated)")
            xxx = pd.DataFrame(fc)
        else:
            xxx = xxx.join(fc, how="outer")

        colors = pd.DataFrame()
        colors[sns.color_palette("colorblind")[0]] = fc.index.str.contains("Wnt")
        colors[sns.color_palette("colorblind")[1]] = fc.index.str.contains("CTRL")
        colors[sns.color_palette("colorblind")[2]] = fc.index.str.contains("Tcr")
        colors[sns.color_palette("colorblind")[3]] = fc.index.str.contains("Ess")
        colors = colors.apply(lambda j: j[j].index.tolist()[0], axis=1).tolist()

        ax.scatter(fc.rank(ascending=False, method="first"), fc, color=colors, alpha=0.5)
        ax.axhline(y=0, color='black', linestyle='--', lw=0.5)
        ax.set_title(re.sub("_stimulated", "", screens[i + 1]))
        ax.set_xlabel("gRNA rank (stimulated / unstimulated)")

    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.diff_condition.svg".format(prefix)), bbox_inches="tight")

    xxx.columns = xxx.columns.str.extract("(.*)_stimulated")
    xxx.to_csv(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.diff_condition.csv".format(prefix)), index=True)


root_dir = "/scratch/lab_bock/shared/projects/crop-seq"
results_dir = os.path.join(root_dir, "results")
sample_annotation = pd.read_csv(os.path.join(root_dir, "metadata/annotation.csv"))

# get guide annotation
guide_annotation = os.path.join(root_dir, "metadata/guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)


# get guide quantification pre screen
counts = ["plasmid_pool_ESS", "plasmid_pool_TCR", "plasmid_pool_WNT"]
pre_screen_counts = pd.DataFrame()
for count in counts:
    df = pd.read_csv(os.path.join("gRNA_counts", count + "_gRNA_count.tsv"), sep="\t")
    df["library"] = count
    pre_screen_counts = pre_screen_counts.append(df)
pre_screen_counts = pd.pivot_table(pre_screen_counts, index="gRNA_name", columns="library", values="count")
pre_screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.pre_screen.csv"), index=True)
pre_screen_counts = pre_screen_counts.apply(lambda x: (x / x.sum()) * 1e4, axis=0)
pre_screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.pre_screen.norm.csv"), index=True)


# get guide quantification mid screen
counts = ["gDNA_HEKclone4", "gDNA_HEKclone6", "gDNA_Jurkat"]
mid_screen_counts = pd.DataFrame()
for count in counts:
    df = pd.read_csv(os.path.join("gRNA_counts", count + "_gRNA_count.tsv"), sep="\t")
    df["library"] = count
    mid_screen_counts = mid_screen_counts.append(df)
mid_screen_counts = pd.pivot_table(mid_screen_counts, index="gRNA_name", columns="library", values="count")
mid_screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.mid_screen.csv"), index=True)
mid_screen_counts = mid_screen_counts.apply(lambda x: (x / x.sum()) * 1e4, axis=0)
mid_screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.mid_screen.norm.csv"), index=True)

# merge output (reads in constructs and assignemnts) of each sample
reads_all = scores_all = assignment_all = pd.DataFrame()

for sample_name in sample_annotation[~sample_annotation["grna_library"].isnull()]["sample_name"].unique():
    reads = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "quantification", "guide_cell_quantification.csv"))
    scores = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "quantification", "guide_cell_scores.csv"), index_col=0).reset_index()
    assignment = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "quantification", "guide_cell_assignment.csv"))

    reads["sample"] = scores["sample"] = assignment["sample"] = sample_name
    reads["experiment"] = scores["experiment"] = assignment["experiment"] = sample_name

    reads_all = reads_all.append(reads)
    scores_all = scores_all.append(scores)
    assignment_all = assignment_all.append(assignment)

reads_all.to_csv(os.path.join(results_dir, "guide_cell_quantification.all.csv"), index=False)
scores_all.to_csv(os.path.join(results_dir, "guide_cell_scores.all.csv"), index=False)
assignment_all.to_csv(os.path.join(results_dir, "guide_cell_assignment.all.csv"), index=False)

reads = pd.read_csv(os.path.join(results_dir, "guide_cell_quantification.all.csv"))
scores = pd.read_csv(os.path.join(results_dir, "guide_cell_scores.all.csv"))
assignment = pd.read_csv(os.path.join(results_dir, "guide_cell_assignment.all.csv"))


screen_counts = pd.pivot_table(assignment.groupby(["experiment", "assignment"]).apply(len).reset_index(), index="assignment", columns="experiment", fill_value=0)
screen_counts.columns = screen_counts.columns.droplevel(level=0)
screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.screen.csv"), index=True)
screen_counts = screen_counts.apply(lambda x: (x / x.sum()) * 1e4, axis=0)
screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.mid_screen.norm.csv"), index=True)


# for grna in mid_screen_counts.index[~mid_screen_counts.index.isin(screen_counts.index)]:
#     screen_counts.loc[grna, ] = 0
pre_sum = (pre_screen_counts.sum(axis=1) / pre_screen_counts.sum().sum()) * 1e4
pre_max = (pre_screen_counts.max(axis=1) / pre_screen_counts.sum().sum()) * 1e4

for s1, s1_ in [
        (pd.DataFrame(pre_screen_counts), "original"),
        (pd.DataFrame(pre_sum), "plasmid_sum"),
        (pd.DataFrame(pre_max), "plasmid_max"),
        (mid_screen_counts, "mid_screen")]:
    for s2, s2_ in [(screen_counts, "crop_screen")]:
        gRNA_scatter(s1, s2, prefix="-".join([s1_, s2_]))
        gRNA_maplot(s1, s2, prefix="-".join([s1_, s2_]))
        gRNA_rank(s1, s2, prefix="-".join([s1_, s2_]))
        gRNA_rank_stimulus(s1, s2, prefix="-".join([s1_, s2_]))
