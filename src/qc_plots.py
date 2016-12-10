#!/usr/bin/env python

import matplotlib
import os
import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from looper.models import Project


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


# Start project, add samples
prj = Project(os.path.join("metadata", "config.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = results_dir = os.path.join("results")


# This is produced from the pipeline statistics report and assembled
# by looper with the command
# `looper summarize metadata/config.yaml`
stats = pd.read_csv("crop-seq_stats_summary.tsv", sep="\t")
gene_thresholds = [500]


# Start gathering
for sample in [s for s in prj.samples if hasattr(s, "replicate")]:  # [s for s in prj.samples if hasattr(s, "replicate")]
    print(sample.name)
    sample_mask = stats['sample_name'] == sample.name
    # Get number of assigned cells
    try:
        sample_assignments = pd.read_csv(
            os.path.join(sample.paths.sample_root, "gRNA_assignment", "guide_cell_assignment.csv")
        ).set_index('cell')
        # total number of assigned cells
        stats.loc[sample_mask, "total_grna_assigned_cells"] = sample_assignments.index.drop_duplicates().shape[0]
        # average number of bps covered per cell
        stats.loc[sample_mask, "mean_grna_basepairs_covered"] = sample_assignments['score'].mean()
        # average number of bps covered per cell
        stats.loc[sample_mask, "mean_concordance_ratio"] = sample_assignments['concordance_ratio'].mean()
    except IOError:
        pass  # if error, it will automatically be pd.np.nan
    try:
        sample_guide_quantification = pd.read_csv(
            os.path.join(sample.paths.sample_root, "gRNA_assignment", "guide_cell_quantification.csv")
        )
        grnas_per_cell = sample_guide_quantification.groupby(['cell'])['molecule'].nunique()
        assignments_per_cell = sample_guide_quantification.groupby(['cell'])['chrom'].nunique()

        # mean grna molecules per cell
        stats.loc[sample_mask, "mean_grna_molecules"] = grnas_per_cell.mean()
        # unambiguous assignemnts
        stats.loc[sample_mask, "percent_unambiguous_grna_assignments"] = ((assignments_per_cell == 1).sum() / float(assignments_per_cell.shape[0])) * 100
    except IOError:
        pass  # if error, it will automatically be pd.np.nan

    for n_genes in gene_thresholds:
        # Gather additional metrics from transcriptome:
        try:
            exp = pd.read_csv(
                os.path.join(sample.paths.sample_root, "digital_expression.{}genes.tsv".format(n_genes)),
                sep="\t").set_index("GENE")
        except IOError:
            continue
        # reads per cell
        reads_per_cell = exp.sum(axis=0)
        stats.loc[sample_mask, "{}genes_total_used_reads".format(n_genes)] = reads_per_cell.sum()

        # # genes per cell
        genes_per_cell = exp.apply(lambda x: (x > 0).sum(), axis=0)
        # # % mitochondrial
        mito_per_cell = (exp.ix[exp.index[exp.index.str.contains("^MT-")]].sum(axis=0) / reads_per_cell) * 100
        # # % ribosomal proteins
        ribo_per_cell = (exp.ix[exp.index[exp.index.str.contains("^RP")]].sum(axis=0) / reads_per_cell) * 100

        # for each, add to stats table distribution values
        for metric in ["reads_per_cell", "genes_per_cell", "mito_per_cell", "ribo_per_cell"]:
            stats.loc[sample_mask, "{}genes_mean_{}".format(n_genes, metric)] = eval(metric).mean()
            stats.loc[sample_mask, "{}genes_median_{}".format(n_genes, metric)] = eval(metric).median()

        # UMI duplication stats
        umi = pd.read_csv(
            os.path.join(sample.paths.sample_root, "cell_umi_barcodes.{}genes.tsv".format(n_genes)),
            sep="\t")

        # get duplication per cell
        dups = umi.groupby(['Cell Barcode'])['Num_Obs'].apply(lambda x: (x == 1).sum())
        duplicates_per_cell = dups / reads_per_cell

        # percent unique UMIs and UMI duplication
        perc_uni_umi = ((umi["Num_Obs"] == 1).sum() / float(umi.shape[0])) * 100
        stats.loc[sample_mask, "{}genes_percent_unique_umis".format(n_genes)] = perc_uni_umi

        # gRNA assignemnt
        # number of cells in transcriptome matrix that are assigned
        assigned = float(sample_assignments.index.isin(exp.columns).sum())
        stats.loc[sample_mask, "{}genes_grna_assigned_cells".format(n_genes)] = assigned
        stats.loc[sample_mask, "{}genes_grna_unassigned_cells".format(n_genes)] = exp.shape[1] - assigned
        stats.loc[sample_mask, "{}genes_percent_grna_assigned_cells".format(n_genes)] = (assigned / exp.shape[1]) * 100

        stats.to_csv(os.path.join(results_dir, "stats.csv"), index=True)

        # Plots

        # plot distributions of transcriptome stats
        fig, axis = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False)
        axis = axis.flatten()
        for i, d in enumerate(["reads_per_cell", "genes_per_cell", "mito_per_cell", "ribo_per_cell"]):
            sns.distplot(eval(d), ax=axis[i])
            axis[i].set_title(d)
        sns.despine(fig)
        fig.savefig(
            os.path.join(results_dir, "stats.{}.{}.distplots.svg".format(sample.name, n_genes)),
            bbox_inches="tight")

        # plot distibution
        fig, axis = plt.subplots(1, figsize=(12, 8))
        sns.distplot(umi["Num_Obs"], kde=False, ax=axis, hist_kws={"log": "True"})
        axis.set_ylabel("Frequency")
        axis.set_xlabel("UMI observations")
        fig.savefig(os.path.join(results_dir, "stats.{}.{}.umi.distplot.svg".format(sample.name, n_genes)), bbox_inches="tight")

        # plot pairwise distributions of each
        combined_stats = pd.concat([
            reads_per_cell, genes_per_cell, mito_per_cell, ribo_per_cell, duplicates_per_cell.ix[reads_per_cell.index]  # transcriptome stats
        ], 1)
        combined_stats.columns = ["reads_per_cell", "genes_per_cell", "mito_per_cell", "ribo_per_cell", "duplicates_per_cell"]
        combined_stats.to_csv(os.path.join(results_dir, "stats.{}.{}genes.all_stats_per_cell.all_cells.csv".format(sample.name, n_genes)))

        g = sns.pairplot(combined_stats, diag_kws={"bins": 100}, plot_kws={"s": 10, "alpha": 0.5})
        g.fig.savefig(
            os.path.join(results_dir, "stats.{}.{}genes.pairwise_dists.svg".format(sample.name, n_genes)),
            bbox_inches="tight")

        # plot pairwise distributions of only assigned cells with good transcriptome
        combined_stats = pd.concat([
            reads_per_cell, genes_per_cell, mito_per_cell, ribo_per_cell, duplicates_per_cell.ix[reads_per_cell.index],  # transcriptome stats
            sample_assignments['score'].ix[reads_per_cell.index], sample_assignments['concordance_ratio'].ix[reads_per_cell.index],  # grna stats for matching cells
            grnas_per_cell.ix[reads_per_cell.index], assignments_per_cell.ix[reads_per_cell.index]  # grna stats for matching cells
        ], 1).dropna()
        combined_stats.columns = [
            "reads_per_cell", "genes_per_cell", "mito_per_cell", "ribo_per_cell", "duplicates_per_cell",
            "grna_bp_covered_per_cell", "grna_concordance_ratio_per_cell", "grna_molecules_per_cell", "detected_grnas_per_cell"]
        combined_stats.to_csv(
            os.path.join(results_dir, "stats.{}.{}genes.all_stats_per_cell.csv".format(sample.name, n_genes)), index=True)

        g = sns.pairplot(combined_stats, diag_kws={"bins": 100}, plot_kws={"s": 10, "alpha": 0.5})
        g.fig.savefig(
            os.path.join(results_dir, "stats.{}.{}genes.pairwise_dists.all_assigned.svg".format(sample.name, n_genes)),
            bbox_inches="tight")
        g.fig.savefig(
            os.path.join(results_dir, "stats.{}.{}genes.pairwise_dists.all_assigned.png".format(sample.name, n_genes)),
            bbox_inches="tight", dpi=300)

        # Number of genes assigned depending on the minimum number of genes required
        for metric in ["genes_per_cell", "reads_per_cell"]:
            total = list()
            assigned = list()
            for j in range(1, 50000, 10):
                pos = eval(metric) > j
                d = eval(metric)[pos]
                try:
                    original_names = pd.Series(d.index.str.split("_")).apply(lambda x: x[0])
                    assigned.append(original_names.isin(sample_assignments.index.dropna()).sum())
                    total.append(pos.sum())
                except AttributeError:
                    break

            fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
            axis[0].plot((range(1, len(total) * 10, 10)), (total), lw=4, color="b", label="All cells")
            axis[0].plot((range(1, len(total) * 10, 10)), (assigned), lw=4, color="r", label="Cells with gRNA assigned")
            axis[1].plot(np.log2(range(1, len(total) * 10, 10)), np.log2(total), lw=4, color="b", label="All cells")
            axis[1].plot(np.log2(range(1, len(total) * 10, 10)), np.log2(assigned), lw=4, color="r", label="Cells with gRNA assigned")
            axis[0].set_xlabel(metric + "")
            axis[1].set_xlabel(metric + " (log2)")
            axis[0].set_ylabel("Cells")
            axis[1].set_ylabel("Cells (log2)")
            axis[1].set_xlim((8, 18))
            plt.legend()
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "stats.{}.{}genes.assigned_cells.varying_{}.svg".format(sample.name, n_genes, metric)), bbox_inches="tight")

            # cumulative sum
            d = pd.DataFrame([total, assigned], index=['total', 'assigned'], columns=range(1, len(total) * 10, 10)).T
            cumsum = d.cumsum() / d.sum()

            fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
            axis[0].plot(cumsum.index, cumsum['total'], lw=4, color="b", label="All cells")
            axis[0].plot(cumsum.index, cumsum['assigned'], lw=4, color="r", label="Cells with gRNA assigned")
            axis[1].plot(np.log2(cumsum.index), cumsum['total'], lw=4, color="b", label="All cells")
            axis[1].plot(np.log2(cumsum.index), cumsum['assigned'], lw=4, color="r", label="Cells with gRNA assigned")
            axis[0].set_xlabel(metric + "")
            axis[1].set_xlabel(metric + " (log2)")
            axis[0].set_ylabel("Cells (fraction of total)")
            axis[1].set_ylabel("Cells (fraction of total)")
            axis[1].set_xlim((8, 18))
            plt.legend()
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "stats.{}.{}.assigned_cells.varying_{}.cumulative.svg".format(sample.name, n_genes, metric)), bbox_inches="tight")
        plt.close('all')


# Plot stats

# transform some metrics to percentages
stats["filtered %"] = 100. - ((stats["FilterBAM in total"] / stats["input_file read1"].astype(float)) * 100)
stats["start_trimmed %"] = 100. - ((stats["TrimStartingSequence in total"] / stats["FilterBAM in total"].astype(float)) * 100)
stats["polyA_trimmed %"] = 100 - ((stats["PolyATrimmer in total"] / stats["TrimStartingSequence in total"].astype(float)) * 100)
stats["bead_fixed_first_base %"] = (stats["FIXED_FIRST_BASE"] / stats["NUM_BEADS"].astype(float)) * 100.
stats["bead_other_error %"] = (stats["OTHER_ERROR_COUNT"] / stats["NUM_BEADS"].astype(float)) * 100.
stats["bead_single_umi_error %"] = (stats["SINGLE_UMI_ERROR"] / stats["NUM_BEADS"].astype(float)) * 100.
stats["bead_primer_match %"] = (stats["PRIMER_MATCH"] / stats["NUM_BEADS"].astype(float)) * 100.
stats["bead_synthesis_error %"] = stats["synthesis_error %"]
stats["DigitalExpression_500genes total_used_reads %"] = (stats["DigitalExpression_500genes total_used_reads"] / stats["input_file read1"].astype(float)) * 100

stats["surviving filtering %"] = (stats["STAR Number of input reads"] / stats["input_file read1"].astype(float)) * 100
stats["Alignment rate"] = stats["STAR Uniquely mapped reads %"] + stats["STAR % of reads mapped to multiple loci"]

odds = False
if 'total_cell_estimation' in stats.columns:
    odds = True
    stats["500genes change_observed_expected number_cells"] = (stats["DigitalExpression_500genes number_cells"] / stats['total_cell_estimation'].astype(float)) * 100


# select some stats

g1 = [
    # initial reads and along triming, filtering
    "input_file read1"]
g2 = [
    "filtered %",
    "start_trimmed %",
    "polyA_trimmed %",
    # mapping
    "STAR % of reads unmapped: too many mismatches",
    "STAR % of reads unmapped: too short",
    "STAR % of reads unmapped: other",
    "STAR % of reads mapped to too many loci",
    "STAR % of reads mapped to multiple loci",
    "STAR Uniquely mapped reads %",

    # bead quality
    "bead_fixed_first_base %",
    "bead_other_error %",
    "bead_single_umi_error %",
    "bead_primer_match %",
    "bead_synthesis_error %",

    # efficiency
    "DigitalExpression_500genes total_used_reads %",
    "500genes percent_unique_umis"]
if odds:
    g2 += ["500genes change_observed_expected number_cells"]

g3 = [
    # transcriptome
    "DigitalExpression_500genes number_cells",
    "DigitalExpression_500genes number_genes",
    "DigitalExpression_500genes reads_per_cell:mean",
    "DigitalExpression_500genes reads_per_cell:median",
    "DigitalExpression_500genes reads_per_cell:std",
    "DigitalExpression_500genes 1reads_to_coverage_:genes_per_cell:mean",
    "DigitalExpression_500genes 1reads_to_coverage_:genes_per_cell:median",
    "DigitalExpression_500genes 1reads_to_coverage_:genes_per_cell:std",

    # grna assignment
    "total_grna_assigned_cells",
    "mean_grna_basepairs_covered",
    "mean_grna_molecules",
    "percent_unambiguous_grna_assignments",
    "500genes_grna_assigned_cells",
    "500genes_grna_unassigned_cells",
    "500genes_percent_grna_assigned_cells"]


metrics = g1 + g2 + g3
stats_sel = stats.T.ix[metrics].astype(float)
stats_sel.to_csv(os.path.join(results_dir, "stats.selected.csv"), index=True)

# plot selected metrics
fig, axis = plt.subplots(3, 1, figsize=(18, 18), gridspec_kw = {'height_ratios': [0.8, 4.1, 4.1]})
sns.heatmap(stats_sel.ix[g1], ax=axis[0], annot=True, fmt='.0f', square=True)
sns.heatmap(stats_sel.ix[g2], ax=axis[1], annot=True, fmt='.1f', square=True)
sns.heatmap(stats_sel.ix[g3], ax=axis[2], annot=True, fmt='.1f', square=True)
for ax in axis[:-1]:
    ax.set_xticklabels(ax.get_xticklabels(), visible=False)
axis[-1].set_xticklabels(axis[-1].get_xticklabels(), rotation=90)
for ax in axis:
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
fig.savefig(os.path.join(results_dir, "stats.selected.heatmap.svg"), bbox_inches="tight")


# only for unmerged samples
sel_samples = [s.name for s in prj.samples if hasattr(s, "replicate") and s.name in stats_sel.columns]

fig, axis = plt.subplots(3, 1, figsize=(18, 18), gridspec_kw = {'height_ratios': [0.8, 4.1, 4.1]})
sns.heatmap(stats_sel[sel_samples].ix[g1], ax=axis[0], annot=True, fmt='.0f', square=True)
sns.heatmap(stats_sel[sel_samples].ix[g2], ax=axis[1], annot=True, fmt='.1f', square=True)
sns.heatmap(stats_sel[sel_samples].ix[g3], ax=axis[2], annot=True, fmt='.1f', square=True)
for ax in axis[:-1]:
    ax.set_xticklabels(ax.get_xticklabels(), visible=False)
axis[-1].set_xticklabels(axis[-1].get_xticklabels(), rotation=90)
for ax in axis:
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
fig.savefig(os.path.join(results_dir, "stats.selected.unmerged.heatmap.svg"), bbox_inches="tight")


fig, axis = plt.subplots(1, figsize=(18, 18))
sns.heatmap(stats_sel[sel_samples].ix[g1 + g2 + g3].apply(
    lambda x: (x - x.min()) / (x.max() - x.min()), axis=1),
    # lambda x: (x - x.mean()) / x.std(), axis=1),
    ax=axis, annot=True, fmt='.2f', square=True)
axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
fig.savefig(os.path.join(results_dir, "stats.selected.unmerged.heatmap.standard_score.svg"), bbox_inches="tight")

#

#

# Figure 1h
# select metrics to show in main figure
metrics = [
    "input_file read1",
    "surviving filtering %",
    "Alignment rate",
    "500genes_percent_unique_umis",
    "DigitalExpression_500genes number_cells",
    "500genes_percent_grna_assigned_cells",
    "DigitalExpression_500genes reads_per_cell:mean",
    "DigitalExpression_500genes 1reads_to_coverage_:genes_per_cell:mean"
]
stats_sel = stats.T.ix[metrics].astype(float)

# get mean of unmerged Jurkat runs
prj.sheet.df.groupby(["experiment", "condition"])

for group, rows in prj.sheet.df.groupby(["experiment", "condition"]):
    means = stats_sel[rows['sample_name'].tolist()].mean(axis=1)
    means["input_file read1"] = stats_sel[rows['sample_name'].tolist()].ix["input_file read1"].sum()
    means["DigitalExpression_500genes number_cells"] = stats_sel[rows['sample_name'].tolist()].ix["DigitalExpression_500genes number_cells"].sum()
    stats_sel["_".join(group) + "_mean"] = means

names = [
    "Read pairs sequenced",
    "Quality filtered",
    "Alignment rate",
    "Unique UMIs (%)",
    u"Cells with ≥500 genes",
    u"Cells assigned (≥500 genes)",
    "Reads per cell (mean)",
    "Genes per cell (mean)",
]
p = stats_sel.ix[metrics].drop_duplicates()
p.index = names

# million
p.ix["Read pairs sequenced"] = p.ix["Read pairs sequenced"] / 1e6

# subset
p = p[p.columns[p.columns.str.contains("HEK|mean")][::-1]]

# plot as bubbles
import matplotlib.patches as mpatches

fig, axis = plt.subplots(1, figsize=(12, 12))
for i, metric in enumerate(p.index):
    x = p.ix[metric] / p.ix[metric].min()
    for j, sample in enumerate(p.columns):
        axis.add_patch(mpatches.Circle((-j, -i), radius=np.sqrt(x[j] / np.pi), ec="none", color=sns.color_palette("colorblind", p.shape[0])[i], alpha=0.5))
        axis.text(-j, (-i - 0.5), s="{0:.1f}".format(p.loc[metric, sample]), horizontalalignment="center")
        if i == 0:
            axis.text(-j, 1, s=sample, horizontalalignment="left", rotation="vertical")
    axis.text(-p.shape[1] - 0.1, -i, s=metric, horizontalalignment="right")
axis.set_xticklabels(axis.get_xticklabels(), visible=False)
axis.set_yticklabels(axis.get_yticklabels(), visible=False)
axis.axis('scaled')
sns.despine(top=True, right=True, left=True, bottom=True)
fig.savefig(os.path.join(results_dir, "figures", "stats.selected.bubbles.svg"), bbox_inches="tight")


#

#


# Figure 1i

# sample_name = "CROP-seq_HEK293T_1_resequenced"
# sample = [s for s in prj.samples if s.name == sample_name][0]

all_counts = pd.DataFrame()
for sample in prj.samples:
    print(sample)
    try:
        exp_counts = pd.read_csv(os.path.join(sample.paths.sample_root, "digital_expression.summary.100genes.tsv"), sep="\t", skiprows=2)
        exp_counts = exp_counts.drop("NUM_TRANSCRIPTS", axis=1).set_index("CELL_BARCODE").squeeze()
        reads = pd.read_csv(os.path.join(sample.paths.sample_root, "gRNA_assignment", "guide_cell_quantification.csv"))
    except IOError:
        continue

    # Assign
    # process
    u = reads.groupby(
        ['cell', 'molecule', 'chrom'])[
        'read_start', 'read_end', 'distance', 'overlap', 'inside', 'mapping_quality', 'strand_agreeement'].max().reset_index()
    # further reduce molecules by solving chromosome conflicts (assign molecule to chromosome with maximum overlap)
    uu = u.ix[u.groupby(['cell', 'molecule']).apply(lambda x: x['overlap'].argmax())]
    # remove no overlaps and reads in wrong strand
    u = uu[uu['strand_agreeement'] == 1]

    # Get a score (sum of bp covered)
    scores = u.groupby(["cell", 'chrom'])['overlap'].sum()
    scores = scores.reset_index().pivot_table("overlap", "cell", "chrom").fillna(0)

    # assign (get max)
    scores["assignment"] = scores.apply(np.argmax, axis=1)
    scores["score"] = scores[list(set(u['chrom']))].apply(max, axis=1)
    # give nan to cells with no overlap (this is because argmax pickus a draw)
    scores.loc[scores['score'] == 0, 'assignment'] = pd.np.nan
    scores.loc[scores['score'] == 0, 'score'] = pd.np.nan

    # Get only assigned cells
    scores = scores.dropna()

    # concordance between reads in same cell
    scores['concordance_ratio'] = scores.drop(["assignment", "score"], axis=1).apply(
        lambda x: x / sum(x), axis=1).max(axis=1)

    counts = pd.DataFrame()
    fig, axis = plt.subplots(3, 3, sharex=True, figsize=(12, 12))
    for i, threshold in enumerate([125, 250, 500, 1000, 2000, 4000, 8000]):
        # Get cells with that many genes
        scores_in = scores[scores.index.isin(exp_counts[exp_counts > threshold].index)]
        print(threshold, scores_in.shape[0])

        s = scores_in.join(exp_counts).dropna()
        axis.flat[i].scatter(np.log2(s['NUM_GENES']), s['concordance_ratio'])
        axis.flat[i].set_title("Cells with at least {} genes.\nMean concordance: {}".format(threshold, s['concordance_ratio'].mean()))

        if scores_in.shape[0] == 0:
            break

        # Quantify
        # in how many cells is the dominant gRNA dominanting by less than 3 times than itself?
        fold = scores_in.drop(["assignment", "score", "concordance_ratio"], axis=1).max(axis=1) / scores_in.drop(["assignment", "score", "concordance_ratio"], axis=1).sum(axis=1)

        counts.loc["not assigned", threshold] = exp_counts[exp_counts > threshold].shape[0] - scores_in.shape[0]
        counts.loc["% not assigned", threshold] = (counts.loc["not assigned", threshold] / float(exp_counts[exp_counts > threshold].shape[0])) * 100
        counts.loc["singlets", threshold] = (fold >= 0.75).sum()
        counts.loc["% singlets", threshold] = ((fold >= 0.75).sum() / float(exp_counts[exp_counts > threshold].shape[0])) * 100
        counts.loc["impurities", threshold] = (fold < 0.75).sum()
        counts.loc["% impurities", threshold] = ((fold < 0.75).sum() / float(exp_counts[exp_counts > threshold].shape[0])) * 100
    counts.columns.name = "genes_covered"
    counts.index.name = "metric"

    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "stats.{}.grna_assignemnt.mean_concordance_vs_genes.various_thresholds.svg".format(sample.name)), bbox_inches="tight")

    fig, axis = plt.subplots(1)
    sns.barplot(
        x="genes_covered",
        y="value",
        hue="metric",
        data=pd.melt(counts[counts.index.str.contains("%")].reset_index(), id_vars=["metric"]),
        ax=axis)
    axis.set_ylim(0, 100)
    axis.set_xlabel("Genes covered per cell")
    axis.set_ylabel("Percentage of all cells")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "figures", "stats.{}.grna_assignemnt.various_thresholds.svg".format(sample.name)), bbox_inches="tight")

    counts["sample_name"] = sample.name
    all_counts = all_counts.append(counts)

# Plot all combined (for reviewer figure 2c)

mean = all_counts[all_counts['sample_name'] != "CROP-seq_HEK293T_1_resequenced"].drop("sample_name", axis=1)

fig, axis = plt.subplots(1)
sns.barplot(
    x="genes_covered",
    y="value",
    hue="metric",
    data=pd.melt(mean[mean.index.str.contains("%")].reset_index(), id_vars=["metric"]),
    ax=axis)
axis.set_ylim(0, 100)
axis.set_xlabel("Genes covered per cell")
axis.set_ylabel("Percentage of all cells")
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "figures", "stats.all_Jurkat_samples.grna_assignemnt.various_thresholds.svg".format(sample.name)), bbox_inches="tight")


# Figure 1i
mean = all_counts.drop("sample_name", axis=1)

fig, axis = plt.subplots(1)
sns.barplot(
    x="genes_covered",
    y="value",
    hue="metric",
    data=pd.melt(mean[mean.index.str.contains("%")].reset_index(), id_vars=["metric"]),
    ax=axis)
axis.set_ylim(0, 100)
axis.set_xlabel("Genes covered per cell")
axis.set_ylabel("Percentage of all cells")
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "figures", "stats.all_samples.grna_assignemnt.various_thresholds.svg".format(sample.name)), bbox_inches="tight")

#

#

# Supplementary Figure 3g

# Read locations within transcripts

# sample_name = "CROP-seq_HEK293T_1_resequenced"
# sample = [s for s in prj.samples if s.name == sample_name][0]

dists = dict()
for sample in [w for w in prj.samples if s.genome == "human"]:

    # read refFlat
    ref = pd.read_csv(os.path.join("spiked_genomes", "hg38_spiked_HEKlibrary", "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat"), sep="\t", header=None)
    ref = ref[range(6)]
    ref.columns = ["gene_id", "transcript_id", "chrom", "strand", "start", "end"]
    ref2 = ref.groupby(["gene_id", "chrom"])["end"].max()

    # bam file
    bam = pysam.AlignmentFile(os.path.join("results_pipeline", sample.name, "star_gene_exon_tagged.clean.bam"))
    # iterate through reads, get read position
    dist = list()
    for i, aln in enumerate(bam):
        if i % 1e5 == 0:
            print(i)
        if (
            (not aln.has_tag("GE")) or
            aln.is_qcfail or  # failed quality (never happens, but for the future)
            aln.is_secondary or
            aln.is_duplicate or
            np.mean(aln.query_alignment_qualities) < 10  # low mapping Q (never happens, but for the future)
        ):
            continue

        # get distance from 3' end
        try:
            dist.append(ref2[aln.get_tag("GE"), aln.reference_name] - aln.reference_end)
        except KeyError:
            continue
    dists[sample.name] = dist
    count = Counter(dist)
    c = pd.Series(count).rolling(window=10).median()

    fig, axis = plt.subplots(1, figsize=(12, 6))
    axis.plot(c.index, c, color="blue")

    axis2 = axis.twinx()
    axis2.plot(c.index, c.cumsum() / c.astype(float).sum(), color="red")

    axis.set_ylabel('Frequency', color='blue')
    axis2.set_ylabel('Cumulative sum', color='red')

    for tl in axis.get_yticklabels():
        tl.set_color('blue')
    for tl in axis2.get_yticklabels():
        tl.set_color('red')

    axis.set_xlim(-10, 10000)
    sns.despine(fig)
    fig.savefig(os.path.join("results", "figures", "stats.{}.read_lengths.svg".format(sample.name)), bbox_inches="tight")


#

#

# Reviewer figure 1b
for i, sample in enumerate([sample for sample in prj.samples if hasattr(sample, "replicate")]):
    s = pd.read_csv(os.path.join(sample.paths.sample_root, "synthesis_statistics.summary.txt"), skiprows=6, sep="\t").set_index("SYNTHESIS_ERROR_BASE")
    s.columns = [sample.name]
    s.loc["Total", sample.name] = float(pd.read_csv(os.path.join(sample.paths.sample_root, "synthesis_statistics.summary.txt"), skiprows=2, sep="\t")["NUM_BEADS"].ix[0])
    s.loc["No Error", sample.name] = float(pd.read_csv(os.path.join(sample.paths.sample_root, "synthesis_statistics.summary.txt"), skiprows=2, sep="\t")["NO_ERROR"].ix[0])

    if i == 0:
        errors = s
    else:
        errors = errors.join(s)

# add macosko
for name in ["Drop-seq_humanmouse_macosko", "Drop-seq_mouse_retina_macosko"]:
    s = pd.read_csv(os.path.join("../", "dropseq_optimizations", "results_pipeline", name, "synthesis_statistics.summary.txt"), skiprows=6, sep="\t").set_index("SYNTHESIS_ERROR_BASE")
    s.columns = [name]
    s.loc["Total", name] = float(pd.read_csv(os.path.join("../", "dropseq_optimizations", "results_pipeline", name, "synthesis_statistics.summary.txt"), skiprows=2, sep="\t")["NUM_BEADS"].ix[0])
    s.loc["No Error", name] = float(pd.read_csv(os.path.join("../", "dropseq_optimizations", "results_pipeline", name, "synthesis_statistics.summary.txt"), skiprows=2, sep="\t")["NO_ERROR"].ix[0])
    errors = errors.join(s)

errors = errors.fillna(0).sort_index()

# transform to percentages
errors_p = errors.apply(lambda x: (x / x[range(1, 9)].sum()) * 100., axis=0)

p = pd.melt(errors_p.reset_index(), id_vars=['SYNTHESIS_ERROR_BASE'])

p.loc[p['SYNTHESIS_ERROR_BASE'].apply(type) == int, 'SYNTHESIS_ERROR_BASE'] = abs(p[p['SYNTHESIS_ERROR_BASE'].apply(type) == int]['SYNTHESIS_ERROR_BASE'] - 9)
p.loc[:, 'SYNTHESIS_ERROR_BASE'] = p['SYNTHESIS_ERROR_BASE'].replace(9, 0)

p.loc[p['variable'].str.contains("_macosko"), "group"] = "Macosko"
p.loc[p['variable'].str.contains("r[1-2]"), "group"] = "Datlinger - old beads"
p.loc[p['variable'].str.contains("r[3-6]"), "group"] = "Datlinger - new beads"

p = p.sort_values(["SYNTHESIS_ERROR_BASE", "group", "variable"], ascending=False)

fig, axis = plt.subplots(1)
sns.barplot(data=p[p['SYNTHESIS_ERROR_BASE'].apply(type) == int], x='SYNTHESIS_ERROR_BASE', y='value', hue='variable', ax=axis)
axis.set_xlabel("Base number")
axis.set_ylabel("Beads with errors in base (%)")
sns.despine(fig)
fig.savefig(os.path.join("results", "figures", "FigR1b.bead_errors.all_sample.svg"), bbox_inches="tight")

fig, axis = plt.subplots(1)
sns.barplot(data=p[p['SYNTHESIS_ERROR_BASE'].apply(type) == int], x='SYNTHESIS_ERROR_BASE', y='value', hue='group', ax=axis)
axis.set_xlabel("Base number")
axis.set_ylabel("Beads with errors in base (%)")
sns.despine(fig)
fig.savefig(os.path.join("results", "figures", "FigR1b.bead_errors.groups.svg"), bbox_inches="tight")

p2 = p.copy()
p2['value'] = 100 - p2['value']

fig, axis = plt.subplots(1)
sns.barplot(data=p2[p2['SYNTHESIS_ERROR_BASE'] == "No Error"], x='SYNTHESIS_ERROR_BASE', y='value', hue='group')
axis.set_xlabel("Experiment group")
axis.set_ylabel("Beads with errors (% of total)")
axis.set_ylim((0, 100))
sns.despine(fig)
fig.savefig(os.path.join("results", "figures", "FigR1b.bead_errors.groups.simple.svg"), bbox_inches="tight")


#

#

# Reviewer figure 1c
reads_per_cell = dict()
genes_per_cell = dict()
for i, sample in enumerate([sample for sample in prj.samples if hasattr(sample, "replicate")]):
    for n_genes in [500]:
        print(sample.name, n_genes)
        # Gather additional metrics from transcriptome:
        try:
            exp = pd.read_csv(
                os.path.join(sample.paths.sample_root, "digital_expression.{}genes.tsv".format(n_genes)),
                sep="\t").set_index("GENE")
        except IOError:
            continue
        # reads per cell
        reads_per_cell[sample.name] = exp.sum(axis=0).values
        # # genes per cell
        genes_per_cell[sample.name] = exp.apply(lambda x: (x > 0).sum(), axis=0).values

# add macosko
for name in ["Drop-seq_humanmouse_macosko", "Drop-seq_mouse_retina_macosko"]:
    for n_genes in [500]:
        print(name, n_genes)
        # Gather additional metrics from transcriptome:
        try:
            exp = pd.read_csv(
                os.path.join("../", "dropseq_optimizations", "results_pipeline", name, "digital_expression.{}genes.tsv".format(n_genes)),
                sep="\t").set_index("GENE")
        except IOError:
            continue
        # reads per cell
        reads_per_cell[name] = exp.sum(axis=0).values
        # # genes per cell
        genes_per_cell[name] = exp.apply(lambda x: (x > 0).sum(), axis=0).values

# add Datlinger uncorrected
old = [
    "CROP-seq_HEK293T_1_resequenced",
    "CROP-seq_Jurkat_TCR_stimulated_r1",
    "CROP-seq_Jurkat_TCR_stimulated_r2",
    "CROP-seq_Jurkat_TCR_unstimulated_r1",
    "CROP-seq_Jurkat_TCR_unstimulated_r2",
    "Drop-seq_HEK293T-3T3"]
for name in old:
    for n_genes in [500]:
        print(name, n_genes)
        # Gather additional metrics from transcriptome:
        try:
            exp = pd.read_csv(
                os.path.join("results_pipeline_uncorrected", name, "digital_expression.{}genes.tsv".format(n_genes)),
                sep="\t").set_index("GENE")
        except IOError:
            continue
        # reads per cell
        reads_per_cell[name + "-uncorrected"] = exp.sum(axis=0).values
        # # genes per cell
        genes_per_cell[name + "-uncorrected"] = exp.apply(lambda x: (x > 0).sum(), axis=0).values


names = pd.Series([q.name for q in prj.samples] + [n + "-uncorrected" for n in old])
groups = {
    "Datlinger - old beads uncorrected": names[names.str.contains("r[1-2]") & names.str.contains("uncorrected")],
    "Datlinger - old beads corrected": names[names.str.contains("r[1-2]") & (~names.str.contains("uncorrected"))]
}

fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 4 * 2))
for group, samples in groups.items():
    sns.distplot(np.concatenate([reads_per_cell[k] for k in samples]), ax=axis[0][0], label=group)
    sns.distplot(np.concatenate([genes_per_cell[k] for k in samples]), ax=axis[0][1], label=group)
    sns.distplot(np.log2(np.concatenate([reads_per_cell[k] for k in samples])), ax=axis[1][0], label=group)
    sns.distplot(np.log2(np.concatenate([genes_per_cell[k] for k in samples])), ax=axis[1][1], label=group)
axis[0][0].set_title("Reads per cell")
axis[0][1].set_title("Genes per cell")
axis[0][0].set_xlabel("Reads per cell")
axis[0][1].set_xlabel("Genes per cell")
axis[1][0].set_xlabel("Reads per cell (log2)")
axis[1][1].set_xlabel("Genes per cell (log2)")
axis[0][0].set_ylabel("Density")
axis[1][0].set_ylabel("Density")
axis[1][1].legend()
sns.despine(fig)
fig.savefig(os.path.join("results", "figures", "FigR1c.dropseq_quality.groups.svg"), bbox_inches="tight")


names = pd.Series([q.name for q in prj.samples])
groups = {
    "Datlinger - old beads": names[names.str.contains("r[1-2]")],
    "Datlinger - new beads": names[names.str.contains("r[3-6]")],
    "Macosko - mix": ["Drop-seq_humanmouse_macosko"],
    "Macosko - retina": ["Drop-seq_mouse_retina_macosko"],
}

fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 4 * 2))
for group, samples in groups.items():
    sns.distplot(np.concatenate([reads_per_cell[k] for k in samples]), ax=axis[0][0], label=group)
    sns.distplot(np.concatenate([genes_per_cell[k] for k in samples]), ax=axis[0][1], label=group)
    sns.distplot(np.log2(np.concatenate([reads_per_cell[k] for k in samples])), ax=axis[1][0], label=group)
    sns.distplot(np.log2(np.concatenate([genes_per_cell[k] for k in samples])), ax=axis[1][1], label=group)
axis[0][0].set_title("Reads per cell")
axis[0][1].set_title("Genes per cell")
axis[0][0].set_xlabel("Reads per cell")
axis[0][1].set_xlabel("Genes per cell")
axis[1][0].set_xlabel("Reads per cell (log2)")
axis[1][1].set_xlabel("Genes per cell (log2)")
axis[0][0].set_ylabel("Density")
axis[1][0].set_ylabel("Density")
axis[1][1].legend()
sns.despine(fig)
fig.savefig(os.path.join("results", "figures", "FigR1c.dropseq_quality.groups.simple.svg"), bbox_inches="tight")


#

#

# Figure SX
stats_sel = pd.read_csv(os.path.join(results_dir, "stats.selected.csv"), index_col=0)

g1 = ["input_file read1"]

g2 = [
    "STAR Uniquely mapped reads %",
    "STAR % of reads mapped to multiple loci",
    "bead_synthesis_error %",
    "500genes_percent_grna_assigned_cells",
    "percent_unambiguous_grna_assignments"]
g3 = [
    "DigitalExpression_500genes number_cells",
    "DigitalExpression_500genes reads_per_cell:mean",
    "DigitalExpression_500genes 1reads_to_coverage_:genes_per_cell:mean",
    "total_grna_assigned_cells",
]

# plot selected metrics
sel_samples = [s.name for s in prj.samples if hasattr(s, "replicate") and s.name in stats_sel.columns]

fig, axis = plt.subplots(3, 1, figsize=(18, 18), gridspec_kw = {'height_ratios': [0.8, 4.1, 4.1]})
sns.heatmap(stats_sel[sel_samples].ix[g1], ax=axis[0], annot=True, fmt='.0f', square=True, cmap="Blues")
sns.heatmap(stats_sel[sel_samples].ix[g2], ax=axis[1], annot=True, fmt='.1f', square=True, cmap="Blues", vmax=100.)
sns.heatmap(stats_sel[sel_samples].ix[g3], ax=axis[2], annot=True, fmt='.1f', square=True, cmap="Blues")
for ax in axis[:-1]:
    ax.set_xticklabels(ax.get_xticklabels(), visible=False)
axis[-1].set_xticklabels(axis[-1].get_xticklabels(), rotation=90)
for ax in axis:
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
fig.savefig(os.path.join(results_dir, "figures", "stats.selected.heatmap.svg"), bbox_inches="tight")
