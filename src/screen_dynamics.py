#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import re
from looper.models import Project


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def filter_gRNAs(df, filter_out=None, prefix=""):
    """
    Filter gRNA table (index) by not included gRNAs
    and gRNAs from opposite library
    """
    filtered = df.copy()

    if filter_out is None:
        filter_out = [  # gRNAs which were not used in the screen
            "CTRL00717", "CTRL00728", "CTRL00783", "CTRL00801",
            "CTRL00851", "CTRL00859", "CTRL00868", "CTRL00872",
            "CTRL00878", "CTRL00881", "CTRL00969", "CTRL00972",
            "CTRL00983"]
        filter_out += [  # gRNAs which were not used in the screen
            "Essential_library_ABL1_1", "Essential_library_ABL1_2", "Essential_library_ABL1_3",
            "Essential_library_MAPK1_1", "Essential_library_MAPK1_2", "Essential_library_MAPK1_3",
            "Essential_library_GATA1_1", "Essential_library_GATA1_2", "Essential_library_GATA1_3",
            "Essential_library_BRCA2_1", "Essential_library_BRCA2_2", "Essential_library_BRCA2_3",
            "Essential_library_PARP1_1", "Essential_library_PARP1_2", "Essential_library_PARP1_3"]

    # Filter non-existing
    df = df.ix[df.index[~df.index.isin(filter_out)]]

    all_ctrls = filtered.index[filtered.index.str.contains("Essential|CTRL")]
    others = all_ctrls[~all_ctrls.isin(filter_out)]
    filtered.loc[others, :] = pd.np.nan

    # Filter opposite library
    for col in df.columns:
        if ('TCR' in col) or ('Jurkat' in col) or ('stimulated' in col) or ('unstimulated' in col):
            df.loc[df.index[df.index.str.contains("Wnt")].tolist(), col] = pd.np.nan
            filtered.loc[filtered.index[filtered.index.str.contains("Tcr")].tolist(), col] = pd.np.nan
        elif ('WNT' in col) or ('HEK' in col):
            df.loc[df.index[df.index.str.contains("Tcr")].tolist(), col] = pd.np.nan
            filtered.loc[filtered.index[filtered.index.str.contains("Wnt")].tolist(), col] = pd.np.nan

    if "plasmid_pool_ESS" in df.columns:
        df = df.drop("plasmid_pool_ESS", axis=1)
        filtered = filtered.drop("plasmid_pool_ESS", axis=1)

    # Vizualize distribution of noise
    fig, axis = plt.subplots(1)
    x = df.values.reshape((np.product(df.shape), 1))
    sns.distplot(np.log2(1 + x[~np.isnan(x)]), ax=axis, label="True gRNAs")
    x = filtered.values.reshape((np.product(filtered.shape), 1))
    try:
        sns.distplot(np.log2(1 + x[~np.isnan(x)]), ax=axis, label="False assignments")
    except ZeroDivisionError:
        print("Screen does not have quantified non-existing gRNAs. Noise cannot be estimated.")
        return df

    # plot 95% percentile of noise
    lower_bound = np.percentile(x[~np.isnan(x)], 95)
    axis.axvline(x=np.log2(1 + lower_bound), color='black', linestyle='--', lw=0.5)
    axis.text(x=np.log2(1 + lower_bound), y=0.5, s="95th percentile = {} cells".format(lower_bound))

    axis.set_xlabel("Number of cells assigned (log2)")
    axis.set_ylabel("Density")
    axis.legend()
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.cell_distribution_noise.{}.svg".format(prefix)), bbox_inches="tight")

    # Filter gRNAs below lower bound
    df = df.where(df > lower_bound, pd.np.nan)

    return df


def normalize_by_total(df):
    """
    Normalize number of gRNA molecules / gRNA-assigned cells by total, experiment-wise.
    """
    return df.apply(lambda x: (x / x.sum()) * 1e4, axis=0)


def screen_zscore(series, axis=None, z_score=False, plot=True):
    """
    Calculate screen z score (difference between positive and negative controls).
    """
    Z = lambda pos, neg: 1 - (3 * (np.std(pos) + np.std(neg)) / (abs(np.mean(pos) - np.mean(neg))))

    if z_score:
        series = (series - series.mean()) / series.std()

    pos = series.ix[series.index[series.index.str.contains("Essential")]]
    neg = series.ix[series.index[series.index.str.contains("CTRL")]]

    z = Z(pos, neg)

    # Plot
    if plot:
        pos.name = None
        neg.name = None
        if axis is None:
            fig, axis = plt.subplots(1)
        sns.distplot(pos, ax=axis, label="positive controls")
        sns.distplot(neg, ax=axis, label="negative controls; screen Z-score = {}".format(z))

    return z


def gRNA_scatter(s1, s2, prefix="", text=False, n_labels=30):
    # Scatter of gRNA change
    fig, axis = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(8, 8))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2)  # .fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        x = x.ix[x.index[~x.index.str.contains("Wnt")]]
        if prefix.startswith("mid_screen-"):
            b = x["gDNA_Jurkat"]
        else:
            b = x["plasmid_pool_TCR"]
        x = x.fillna(0)
        b = b.fillna(0)

        colors = pd.DataFrame()
        colors[sns.color_palette("colorblind")[0]] = x.index.str.contains("Wnt")
        colors[sns.color_palette("colorblind")[1]] = x.index.str.contains("CTRL")
        colors[sns.color_palette("colorblind")[2]] = x.index.str.contains("Tcr")
        colors[sns.color_palette("colorblind")[3]] = x.index.str.contains("Ess")
        colors = colors.apply(lambda x: x[x].index.tolist()[0], axis=1).tolist()

        axis[i].scatter(np.log2(1 + x[screen]), np.log2(1 + b), color=colors, alpha=0.5)
        if text:
            for j in x[x.index.str.contains("ETS1|GATA3|RUNX1")].index:
                axis[i].text(np.log2(1 + x[screen].ix[j]), np.log2(1 + b.ix[j]), j)

        # x = y line
        lims = [np.nanmin([np.log2(1 + x[screen]), np.log2(1 + b)]), np.nanmax([np.log2(1 + x[screen]), np.log2(1 + b)])]
        axis[i].plot((lims[0], lims[1]), (lims[0], lims[1]), linestyle='--', color='black', alpha=0.75)

        axis[i].set_title(screen)
    for i in range(0, len(axis), 2):
        axis[i].set_ylabel("gRNA frequency in plasmid (log2)")
    for ax in axis[-2:]:
        ax.set_xlabel("gRNA frequency in CROP-seq screen (log2)")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.scatter.{}svg".format(prefix, "text." if text else "")), bbox_inches="tight")
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.scatter.{}pdf".format(prefix, "text." if text else "")), bbox_inches="tight")


def gRNA_maplot(s1, s2, prefix="", text=False, n_labels=30):
    # Rank of gRNA change
    fig, axis = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8, 8))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2)  # .fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        x = x.ix[x.index[~x.index.str.contains("Wnt")]]
        if prefix.startswith("mid_screen-"):
            b = x["gDNA_Jurkat"]
        else:
            b = x["plasmid_pool_TCR"]
        x = x.fillna(0)
        b = b.fillna(0)

        M = np.log2(x[screen] * b) / 2.
        M = M.replace({-np.inf: 0, np.inf: 9})
        fc = np.log2(1 + x[screen]) - np.log2(1 + b)

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
        if text:
            for j in x[x.index.str.contains("ETS1|GATA3|RUNX1")].index:
                axis[i].text(
                    M.ix[j],
                    fc.ix[j],
                    j)

        axis[i].axhline(y=0, color='black', linestyle='--', lw=0.5)

        axis[i].set_title(screen)

    for i in range(0, len(axis), 2):
        axis[i].set_ylabel("M")
    for ax in axis[-2:]:
        ax.set_xlabel("A")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.maplot.{}svg".format(prefix, "text." if text else "")), bbox_inches="tight")
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.maplot.{}pdf".format(prefix, "text." if text else "")), bbox_inches="tight")


def gRNA_rank(s1, s2, prefix="", text=False, n_labels=30):
    # Rank of gRNA change
    fig, axis = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8, 8))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2)  # .fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        x = x.ix[x.index[~x.index.str.contains("Wnt")]]
        if prefix.startswith("mid_screen-"):
            b = x["gDNA_Jurkat"]
        else:
            b = x["plasmid_pool_TCR"]

        x = x.fillna(0)
        b = b.fillna(0)

        fc = np.log2(1 + x[screen]) - np.log2(1 + b)

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
        if text:
            for j in x[x.index.str.contains("ETS1|GATA3|RUNX1")].index:
                axis[i].text(
                    fc.rank(ascending=False, method="first").ix[j],
                    fc.ix[j],
                    j)
        axis[i].axhline(y=0, color='black', linestyle='--', lw=0.5)

        axis[i].set_title(screen)

    for i in range(0, len(axis), 2):
        axis[i].set_ylabel("gRNA fold-change")
    for ax in axis[-2:]:
        ax.set_xlabel("gRNA rank")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.{}svg".format(prefix, "text." if text else "")), bbox_inches="tight")
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.{}pdf".format(prefix, "text." if text else "")), bbox_inches="tight")

    # Save ranked list
    xx.to_csv(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.csv".format(prefix)), index=True)

    # Save ranked list of gene-level measurements, reduced by mean and min
    m = pd.merge(xx.reset_index(), guide_annotation[["oligo_name", "gene"]], left_on="gRNA_name", right_on="oligo_name").drop("oligo_name", axis=1).set_index(["gene", "gRNA_name"])
    m.groupby(level=[0]).mean().to_csv(os.path.join(results_dir, "gRNA_counts.norm.{}.gene_mean.rank.csv".format(prefix)), index=True)
    m.groupby(level=[0]).min().to_csv(os.path.join(results_dir, "gRNA_counts.norm.{}.gene_min.rank.csv".format(prefix)), index=True)


def gRNA_rank_stimulus(xx, s2, prefix=""):
    # Difference between unstimulated/stimulated
    fig, axis = plt.subplots(1, 3, sharex=False, sharey=True, figsize=(12, 3))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        x = s1.join(s2)  # .fillna(0)
        x = x.iloc[np.random.permutation(len(x))]

        if ("TCR" in screen) or ("Jurkat" in screen):
            x = x.ix[x.index[~x.index.str.contains("Wnt")]]
            if prefix.startswith("mid_screen-"):
                b = x["gDNA_Jurkat"]
            else:
                b = x["plasmid_pool_TCR"]
        elif ("WNT" in screen) or ("HEK" in screen):
            x = x.ix[x.index[~x.index.str.contains("Tcr")]]
            if prefix.startswith("mid_screen-"):
                if "_4_" in prefix:
                    b = x["gDNA_HEKclone4"]
                else:
                    b = x["gDNA_HEKclone6"]
            else:
                b = x["plasmid_pool_WNT"]
        fc = np.log2(1 + x[screen]) - np.log2(1 + b)

        fc.name = screen
        if i == 0:
            xx = pd.DataFrame(fc)
        else:
            xx = xx.join(fc, how="outer")

    screens = s2.columns[::-1]
    for i in range(0, len(s2.columns), 2):
        fc = (xx[screens[i + 1]] - xx[screens[i]]).dropna()

        fc.name = screens[i + 1]
        if i == 0:
            axis[i].set_ylabel("gRNA fold-change (stimulated / unstimulated)")
            xxx = pd.DataFrame(fc)
        else:
            xxx = xxx.join(fc, how="outer")

        colors = pd.DataFrame()
        colors[sns.color_palette("colorblind")[0]] = fc.index.str.contains("Wnt")
        colors[sns.color_palette("colorblind")[1]] = fc.index.str.contains("CTRL")
        colors[sns.color_palette("colorblind")[2]] = fc.index.str.contains("Tcr")
        colors[sns.color_palette("colorblind")[3]] = fc.index.str.contains("Ess")
        colors = colors.apply(lambda j: j[j].index.tolist()[0], axis=1).tolist()

        axis[i].scatter(fc.rank(ascending=False, method="first"), fc, color=colors, alpha=0.5)
        axis[i].axhline(y=0, color='black', linestyle='--', lw=0.5)
        axis[i].set_title(re.sub("_stimulated", "", screens[i + 1]))
        axis[i].set_xlabel("gRNA rank (stimulated / unstimulated)")

    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.diff_condition.svg".format(prefix)), bbox_inches="tight")

    xxx.columns = xxx.columns.str.extract("(.*)_stimulated")
    xxx.to_csv(os.path.join(results_dir, "gRNA_counts.norm.{}.rank.diff_condition.csv".format(prefix)), index=True)


def gRNA_swarmplot(s1, s2, prefix=""):
    # Rank of gRNA change
    fig, axis = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8, 8))
    axis = axis.flatten()

    for i, screen in enumerate(s2.columns[::-1]):
        s = s1.join(s2)  # .fillna(0)
        s = s.iloc[np.random.permutation(len(s))]

        if ("TCR" in screen) or ("Jurkat" in screen) or ("stimulated" in screen) or ("unstimulated" in screen):
            s = s.ix[s.index[~s.index.str.contains("Wnt")]]
            if prefix.startswith("mid_screen-"):
                b = s["gDNA_Jurkat"]
            else:
                b = s["plasmid_pool_TCR"]
            x = s.ix[s.index[s.index.str.contains("Tcr")]]
            y = s.ix[s.index[s.index.str.contains("Essential")]]
            z = s.ix[s.index[s.index.str.contains("CTRL")]]
            b_x = b.ix[s.index[s.index.str.contains("Tcr")]]
            b_y = b.ix[s.index[s.index.str.contains("Essential")]]
            b_z = b.ix[s.index[s.index.str.contains("CTRL")]]
        elif ("WNT" in screen) or ("HEK" in screen):
            s = s.ix[s.index[~s.index.str.contains("Tcr")]]
            if prefix.startswith("mid_screen-"):
                if "_4_" in prefix:
                    b = s["gDNA_HEKclone4"]
                else:
                    b = s["gDNA_HEKclone6"]
            else:
                b = s["plasmid_pool_WNT"]
            x = s.ix[s.index[s.index.str.contains("Wnt")]]
            y = s.ix[s.index[s.index.str.contains("Essential")]]
            z = s.ix[s.index[s.index.str.contains("CTRL")]]
            b_x = b.ix[s.index[s.index.str.contains("Wnt")]]
            b_y = b.ix[s.index[s.index.str.contains("Essential")]]
            b_z = b.ix[s.index[s.index.str.contains("CTRL")]]

        fc_x = np.log2(1 + x[screen]) - np.log2(1 + b_x)
        fc_y = np.log2(1 + y[screen]) - np.log2(1 + b_y)
        fc_z = np.log2(1 + z[screen]) - np.log2(1 + b_z)

        fc_x.name = screen
        fc_y.name = "Essential"
        fc_z.name = "CTRL"

        sns.violinplot(x="variable", y="value", alpha=0.1, inner="box", data=pd.melt(pd.DataFrame([fc_x, fc_y, fc_z]).T), ax=axis[i])
        sns.swarmplot(x="variable", y="value", alpha=0.5, data=pd.melt(pd.DataFrame([fc_x, fc_y, fc_z]).T), ax=axis[i])
        axis[i].axhline(y=0, color='black', linestyle='--', lw=0.5)

        axis[i].set_title(screen)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.norm.{}.violin_swarmplot.svg".format(prefix)), bbox_inches="tight")


prj = Project(os.path.join("metadata", "config.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = results_dir = os.path.join("results", "screen_dynamics")

sample_annotation = prj.sheet.df

# get guide annotation
guide_annotation = os.path.join("metadata", "guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)


experiment = "CROP-seq_Jurkat_TCR"


# get guide quantification pre screen
counts = ["plasmid_pool_ESS", "plasmid_pool_TCR", "plasmid_pool_WNT"]
pre_screen_counts = pd.DataFrame()
for count in counts:
    df = pd.read_csv(os.path.join("gRNA_counts", count + "_gRNA_count.tsv"), sep="\t")
    df["library"] = count
    pre_screen_counts = pre_screen_counts.append(df)
pre_screen_counts = pd.pivot_table(pre_screen_counts, index="gRNA_name", columns="library", values="count")
pre_screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.pre_screen.csv"), index=True)

# get guide quantification mid screen
counts = ["gDNA_HEKclone4", "gDNA_HEKclone6", "gDNA_Jurkat"]
mid_screen_counts = pd.DataFrame()
for count in counts:
    df = pd.read_csv(os.path.join("gRNA_counts", count + "_gRNA_count.tsv"), sep="\t")
    df["library"] = count
    mid_screen_counts = mid_screen_counts.append(df)
mid_screen_counts = pd.pivot_table(mid_screen_counts, index="gRNA_name", columns="library", values="count")
mid_screen_counts.to_csv(os.path.join(results_dir, "gRNA_counts.mid_screen.csv"), index=True)


# get guide quantification from CROP-seq
# merge output (reads in constructs and assignemnts) of each sample

reads = pd.read_csv(os.path.join("results", "{}.guide_cell_gRNA_assignment.all.csv".format(experiment)))
scores = pd.read_csv(os.path.join("results", "{}.guide_cell_scores.all.csv".format(experiment)))
assignment = pd.read_csv(os.path.join("results", "{}.guide_cell_assignment.all.csv".format(experiment)))

screen_counts = pd.pivot_table(assignment.groupby(["experiment", "condition", "assignment"]).apply(len).reset_index(), index="assignment", columns="condition", fill_value=0)
screen_counts.columns = screen_counts.columns.droplevel(level=0)
screen_counts.to_csv(os.path.join(results_dir, "{}.gRNA_counts.screen.csv".format(experiment)), index=True)


# Filter out non-existing gRNAs and opposite gRNA library.
# Use distribution of cells assigned with non-existing gRNAs
# to get 95th percentile of distribution as lower bound of sensitivity
# and further use it to filter out gRNAs with as many cells
pre_screen_counts = filter_gRNAs(pre_screen_counts, prefix="plasmid")
mid_screen_counts = filter_gRNAs(mid_screen_counts, prefix="mid_screen")
screen_counts = filter_gRNAs(screen_counts, prefix="screen")

# Normalize
pre_screen_counts = pre_screen_counts.apply(lambda x: (x / x.sum(skipna=True)) * 1e4, axis=0)
pre_screen_counts.to_csv(os.path.join(results_dir, "{}.gRNA_counts.pre_screen.norm.csv".format(experiment)), index=True)

mid_screen_counts = mid_screen_counts.apply(lambda x: (x / x.sum(skipna=True)) * 1e4, axis=0)
mid_screen_counts.to_csv(os.path.join(results_dir, "{}.gRNA_counts.mid_screen.norm.csv".format(experiment)), index=True)

screen_counts = screen_counts.apply(lambda x: (x / x.sum(skipna=True)) * 1e4, axis=0)
screen_counts.to_csv(os.path.join(results_dir, "{}.gRNA_counts.screen.norm.csv".format(experiment)), index=True)


pre_screen_counts = pd.read_csv(os.path.join(results_dir, "{}.gRNA_counts.pre_screen.norm.csv".format(experiment)), index_col=0)
mid_screen_counts = pd.read_csv(os.path.join(results_dir, "{}.gRNA_counts.mid_screen.norm.csv".format(experiment)), index_col=0)
screen_counts = pd.read_csv(os.path.join(results_dir, "{}.gRNA_counts.screen.norm.csv".format(experiment)), index_col=0)

for s1, s1_ in [
        (pre_screen_counts, "original"),
        (mid_screen_counts, "mid_screen")]:
    for s2, s2_ in [
            (mid_screen_counts, "mid_screen"),
            (screen_counts, "crop_screen")]:
        if s1_ == s2_:
            continue
        # gRNA_rank(s1, s2, prefix="-".join([s1_, s2_]))
        # gRNA_scatter(s1, s2, prefix="-".join([s1_, s2_]))
        # gRNA_maplot(s1, s2, prefix="-".join([s1_, s2_]))
        # # gRNA_rank_stimulus(s1, s2, prefix="-".join([s1_, s2_]))
        gRNA_swarmplot(s1, s2, prefix="-".join([s1_, s2_]))
        # with text labels
        gRNA_rank(s1, s2, prefix="-".join([s1_, s2_]), text=True)
        gRNA_scatter(s1, s2, prefix="-".join([s1_, s2_]), text=True)
        gRNA_maplot(s1, s2, prefix="-".join([s1_, s2_]), text=True)


# Get screen sensitivity
# Z-score
zez = list()

for screen, name in [
        (pre_screen_counts, "original"),
        (mid_screen_counts, "mid_screen"),
        (screen_counts, "crop_screen")]:

    fig, axis = plt.subplots(
        2 if name == "original" else 3, 2 if name == "crop_screen" else 1,
        sharex=False, sharey=False, figsize=(8 if name == "original" else 12, 12 if name == "crop_screen" else 8,))
    axis = axis.flatten()

    for i, col in enumerate(screen.columns[::-1]):
        z = screen_zscore(screen[col].dropna(), axis=axis[i])

        zez.append([name, col, z])

        axis[i].set_title(" ".join([name, col]))
        axis[i].set_xlabel("Number of cells assigned (log2)")
        axis[i].set_ylabel("Density")
        axis[i].legend()
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "gRNA_counts.screen_sensitivity.{}.svg".format(name)), bbox_inches="tight")

zez = pd.DataFrame(zez)
zez.columns = ["timepoint", "sample", "z_score"]
zez["id"] = zez["timepoint"] + " " + zez["sample"]
zez["efficiency"] = 1 / -zez["z_score"]
zez.to_csv(os.path.join(results_dir, "gRNA_counts.screen_sensitivity.z_score.csv"), index=False)

# Barplot across experiments
zez = zez.sort_values(["efficiency"], ascending=False)
fig, axis = plt.subplots(1)
sns.barplot(zez["efficiency"], zez["id"], orient="h", ax=axis)
axis.set_title("Screen sensitivity")
axis.set_xlabel("Sensitivity (1 / Z score)")
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "gRNA_counts.screen_sensitivity.barplot.svg"), bbox_inches="tight")


# plot replicates

# Library composition
counts = pre_screen_counts.join(mid_screen_counts).join(screen_counts)
counts = counts.ix[counts.index[~counts.index.str.contains("Wnt")]].sort_index(ascending=False)
counts = counts[counts.columns[~counts.columns.str.contains("WNT|HEK")]]

# heatmap
fig, axis = plt.subplots(1)
sns.heatmap(
    counts,
    cmap="BuGn",
    ax=axis)
axis.set_xticklabels(axis.get_xticklabels(), rotation="vertical")
axis.set_yticklabels(axis.get_yticklabels(), rotation="horizontal")
fig.savefig(os.path.join(results_dir, "gRNA_counts.fold_distribution.heatmap.svg"), bbox_inches="tight")

# heatmap + z_score
fig, axis = plt.subplots(1)
sns.heatmap(
    counts.apply(lambda x: (x - x.mean()) / x.std(), axis=1),
    ax=axis)
axis.set_xticklabels(axis.get_xticklabels(), rotation="vertical")
axis.set_yticklabels(axis.get_yticklabels(), rotation="horizontal")
fig.savefig(os.path.join(results_dir, "gRNA_counts.fold_distribution.heatmap.z_score.svg"), bbox_inches="tight")

# ranks
fig, axis = plt.subplots(2, 2, sharex=False, sharey=True)
axis = axis.flatten()
for i, name in enumerate(counts.columns):
    axis[i].plot(
        counts[name].sort_values(), (np.arange(len(counts[name])) / float(len(counts[name]))) * 100
    )
    p = np.percentile(counts[name].dropna().sort_values(), 90), np.percentile(counts[name].dropna().sort_values(), 10)
    axis[i].axvline(p[1], linestyle="--", color="black")
    axis[i].axvline(p[0], linestyle="--", color="black")
    axis[i].set_title(name)
    axis[i].text(0.5, 0.5, str(p[0] / p[1]))
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "gRNA_counts.fold_distribution.lineplot.svg"), bbox_inches="tight")


# Plot gRNA abundance variability
# get accuracy prediction
counts_annot = pd.merge(counts.reset_index(), guide_annotation, left_on="gRNA_name", right_on="oligo_name")

# group grnas by gene
# compute amplitude/var/std
var = counts_annot.groupby("gene")[counts.columns.tolist() + ["efficiency_score"]].std()
mean = counts_annot.groupby("gene")[counts.columns.tolist() + ["efficiency_score"]].mean()
varmean = var / mean

# plot power vs variability
fig, axis = plt.subplots(2, 2, sharex=False, sharey=True)
for i, s in enumerate(counts.columns):
    axis.flat[i].scatter(mean[s], var[s])
    axis.flat[i].set_title(s)
    axis.flat[i].set_xlabel("Mean")
    axis.flat[i].set_ylabel("Std")
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "gRNA_counts.mean_std.scatter.svg"), bbox_inches="tight")

# plot predicted variability vs observed
fig, axis = plt.subplots(2, 2, sharex=False, sharey=True)
for i, s in enumerate(counts.columns):
    axis.flat[i].scatter(varmean[s], var["efficiency_score"])
    axis.flat[i].set_title(s)
    axis.flat[i].set_xlabel("qv (Cell abundance)")
    axis.flat[i].set_ylabel("qv (gRNA efficiency)")
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "gRNA_counts.qv-crop_efficiency.scatter.svg"), bbox_inches="tight")
