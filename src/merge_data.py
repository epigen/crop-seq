#!/usr/bin/env python

import argparse
import os
import pandas as pd
import pysam
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def plot_reads_in_constructs(reads):
    from collections import Counter
    # Inspect
    fig, axis = plt.subplots(2, sharex=True)
    # number of barcode reads per cell
    sns.distplot(np.log2(1 + reads.groupby(["cell"]).apply(len)), ax=axis[0], kde=False)
    axis[0].set_xlabel("Reads (log2(1 + x))")
    # number of unique barcode reads per cell
    sns.distplot(np.log2(1 + reads.groupby(["cell"])['molecule'].apply(set).apply(len)), ax=axis[1], kde=False)
    axis[1].set_xlabel("Molecules (log2(1 + x))")
    fig.savefig(os.path.join(results_dir, "barcodes_per_cell.svg"), bbox_inches="tight")

    # index
    u = reads.ix[reads[["chrom", "cell", 'molecule']].drop_duplicates().index]

    # number of unique barcode reads saying inside per cell
    fig, axis = plt.subplots(1)
    sns.distplot(np.log2(1 + u.groupby(["cell"])['molecule'].apply(len)), ax=axis, kde=False)
    axis.set_xlabel("Molecules (log2(1 + x))")
    fig.savefig(os.path.join(results_dir, "barcodes_per_cell.inside.svg"), bbox_inches="tight")

    # concordance between reads in same cell
    inside_ratio = u.groupby(["cell"]).apply(lambda x: sum(x['inside'] == 1) / float(len(x)))
    fig, axis = plt.subplots(1)
    sns.distplot(inside_ratio, kde=False)
    axis.set_xlabel("Ratio molecules overlap gRNA / total")
    fig.savefig(os.path.join(results_dir, "barcodes_per_cell.concordance.svg"), bbox_inches="tight")

    # how many cells are fully out?
    inside_ratio[inside_ratio == 0.0].shape[0]
    # inside?
    inside = inside_ratio[1.0 == inside_ratio]
    inside.shape[0]
    # in-between?
    between = inside_ratio[(1.0 > inside_ratio) & (inside_ratio > 0.0)]
    between.shape[0]

    # distribution of guides for cells inside
    fig, axis = plt.subplots(1)
    c = pd.Series(Counter(u[u['cell'].isin(inside.reset_index()['cell'])][['chrom', 'cell']].drop_duplicates()['chrom'])).sort_values()
    sns.barplot(c.index, c.values, ax=axis)
    axis.set_title("Cells with only molecules inside")
    fig.savefig(os.path.join("barcodes_per_cell.guides_inside.svg"), bbox_inches="tight")

    # distribution of reads regarding constructs (for each guide)
    g = sns.FacetGrid(u, col="chrom", sharex=False, sharey=False)
    g.map(sns.distplot, 'distance', kde=False)
    g.fig.savefig(os.path.join(results_dir, "barcodes_per_cell.distance.svg"), bbox_inches="tight")

    # For all guides plot distribution inside and out
    fig, axis = plt.subplots(2, len(set(u["chrom"])), sharex=False, sharey=False, figsize=(16, 8))
    axis = iter(axis.flatten())
    for inside in [1.0, 0.0]:
        for chrom in set(u["chrom"]):
            ax = axis.next()
            ax.set_title(chrom)
            ax.set_ylabel("Inside" if inside else "Outside")
            subset = u[(u["chrom"] == chrom) & (u["inside"] == inside)]
            if subset.shape[0] > 1:
                sns.distplot(subset['distance'], kde=False, ax=ax)
    fig.savefig(os.path.join(results_dir, "barcodes_per_cell.distance.in_out.svg"), bbox_inches="tight")


def make_assignment(reads, guide_annotation):
    # Assign
    # unique reads per cell
    u = reads.ix[reads[["cell", 'molecule']].drop_duplicates().index]
    # filter out reads in wrong strand
    u = u[u['strand_agreeement'] == 1]
    # get unique reads per cell
    u = reads.ix[reads[["cell", 'molecule']].drop_duplicates().index]

    # Get a score (sum of bp covered)
    scores = u.groupby(["cell", 'chrom'])['overlap'].sum()
    scores = scores.reset_index().pivot_table("overlap", index="cell", columns="chrom").fillna(0)

    # assign (get max)
    scores["assignment"] = scores.apply(np.argmax, axis=1)
    scores["score"] = scores[scores.columns.drop("assignment")].max(axis=1)
    # give nan to cells with no overlap (this is because argmax pickus a draw)
    scores.loc[scores['score'] == 0, 'assignment'] = pd.np.nan
    scores.loc[scores['score'] == 0, 'score'] = pd.np.nan

    # Get only assigned cells
    scores = scores.dropna()

    # Convert to coverage in X times (divide by length of gRNA)
    coverage = scores.drop(['assignment', 'score'], axis=1).apply(lambda x: x / float(len(guide_annotation[guide_annotation["oligo_name"] == x.name]["sequence"].squeeze())), axis=0)
    coverage["maxscore"] = coverage.apply(max, axis=1)
    coverage["assignment"] = coverage.drop("maxscore", axis=1).apply(np.argmax, axis=1)

    cov = coverage[~(coverage["maxscore"] == 0)][["assignment", "maxscore"]]
    cov["maxscore"] = np.log2(1 + cov["maxscore"])

    # Plot bp or covereage covered per cell
    fig, axis = plt.subplots(2, sharex=False, sharey=False)
    sns.distplot(np.log2(scores['score']), kde=False, ax=axis[0])
    axis[0].set_xlabel("gRNA coverage (bp)")
    axis[0].set_ylabel("Frequency")

    sns.distplot(np.log2(cov['maxscore'] + 1), kde=False, ax=axis[1])
    axis[1].set_xlabel("gRNA coverage (X)")
    axis[1].set_ylabel("Frequency")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "barcodes_per_cell.coverage.svg"), bbox_inches="tight")

    # concordance between reads in same cell
    concordance_ratio = (scores["score"] / (scores.drop(["score", "assignment"], axis=1).sum(axis=1))).replace({0.0: pd.np.nan})
    percent = ((concordance_ratio != 1.).sum() / float(concordance_ratio.shape[0])) * 100

    fig, axis = plt.subplots(1, sharex=False, sharey=False)
    sns.distplot(concordance_ratio, kde=False, ax=axis)
    axis.text(0.5, 1000, "{}% non-absolute concordance".format(percent), ha='center', va='center')
    axis.set_xlabel("gRNA coverage (bp)")
    axis.set_ylabel("Frequency")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "barcodes_per_cell.score_concordance_ratio.svg"), bbox_inches="tight")

    # Get assigned cells
    assignment["assigned_gene"] = assignment["assignment"].str.extract("(.*)_\d$")
    # plot assignment stats
    counts = assignment.groupby(['experiment', 'assignment']).apply(len)

    for library in counts.index.levels[0]:
        c = counts.ix[library].sort_values(ascending=False)

        if "TCR" in library:
            c = c.ix[c.index[~c.index.str.contains("Wnt")]]
        else:
            c = c.ix[c.index[~c.index.str.contains("Tcr")]]

        fig, axis = plt.subplots(1, figsize=(6, 12))
        sns.barplot(c.values, c.index, orient="h", ax=axis)
        axis.set_xlabel("gRNA frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "cells_per_barcode.gRNA_assigned.{}.svg".format(library)), bbox_inches="tight")

    counts = assignment.groupby(['experiment', 'assigned_gene']).apply(len)

    for library in counts.index.levels[0]:
        c = counts.ix[library].sort_values(ascending=False)

        if "TCR" in library:
            c = c.ix[c.index[~c.index.str.contains("Wnt")]]
        else:
            c = c.ix[c.index[~c.index.str.contains("Tcr")]]

        fig, axis = plt.subplots(1, figsize=(6, 12))
        sns.barplot(c.values, c.index, orient="h", ax=axis)
        axis.set_xlabel("gRNA frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "cells_per_barcode.gene_assigned.{}.svg".format(library)), bbox_inches="tight")

    return scores, assignment


def big_heatmap(x, assignment):
    from scipy.spatial import distance
    from scipy.cluster import hierarchy
    import sys
    sys.setrecursionlimit(10000)

    # Compute but don't plot first dendrogram.
    col_pairwise_dists = distance.pdist(x.T, metric="euclidean")
    col_linkage = hierarchy.linkage(col_pairwise_dists, method="complete")
    col_dendrogram = hierarchy.dendrogram(col_linkage, no_plot=True)
    # Compute but don't plot second dendrogram.
    row_pairwise_dists = distance.pdist(x, metric="euclidean")
    row_linkage = hierarchy.linkage(row_pairwise_dists, method="single")
    row_dendrogram = hierarchy.dendrogram(row_linkage, no_plot=True)
    # Reorder matrix rows and cols
    x_plot = x.iloc[row_dendrogram['leaves'], col_dendrogram['leaves']]

    # map gRNA to integer/color
    cell_grna_map = assignment.set_index("cell")['assignment'].to_dict()
    integer_map = dict([(val, i) for i, val in enumerate(set(assignment['assignment'].tolist() + ["unnassigned"]))])
    colors = pd.Series(x_plot.columns).apply(lambda x: integer_map[cell_grna_map[x]] if x in cell_grna_map.keys() else integer_map["unnassigned"])

    cell_grna_score = assignment.set_index("cell")['score'].to_dict()
    scores = pd.Series(x_plot.columns).apply(lambda x: cell_grna_score[x] if x in cell_grna_score.keys() else pd.np.nan)

    # plot
    fig = plt.figure(figsize=(8, 8))

    # Plot distance matrix.
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    im1 = ax1.imshow(x_plot, interpolation='nearest', aspect='auto', origin='lower', vmin=x.min().min(), vmax=x.max().max(), cmap="Blues")
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_ylabel("%i genes" % x_plot.shape[0])
    ax1.set_xlabel("%i cells" % x_plot.shape[1])

    # plot gRNA assignment
    ax2 = fig.add_axes([0.1, 0.91, 0.8, 0.015])
    im2 = ax2.imshow(np.array([colors.values]), interpolation='nearest', aspect='auto', cmap="viridis")
    ax2.set_xticks([])
    ax2.set_yticks([])

    # plot gRNA assignment score
    ax3 = fig.add_axes([0.1, 0.925, 0.8, 0.015])
    im3 = ax3.imshow(np.log2(np.array([scores.values])), interpolation='nearest', aspect='auto')
    ax3.set_xticks([])
    ax3.set_yticks([])

    ax4 = fig.add_axes([0.1, 0.05, 0.2, 0.2])
    plt.colorbar(im3, ax=ax4, orientation="horizontal", label='gRNA assignment score')
    ax4.set_axis_off()

    ax5 = fig.add_axes([0.7, 0.05, 0.2, 0.2])
    plt.colorbar(im2, ax=ax5, orientation="horizontal", label='gRNA')
    ax5.set_axis_off()
    im2.colorbar.set_ticks(integer_map.values())
    im2.colorbar.set_ticklabels(integer_map.keys())

    ax6 = fig.add_axes([0.4, 0.05, 0.2, 0.2])
    plt.colorbar(im1, ax=ax6, orientation="horizontal", label='log2(expression)')
    ax6.set_axis_off()

    return fig


root_dir = "/scratch/lab_bock/shared/projects/crop-seq"
results_dir = os.path.join(root_dir, "results")
sample_annotation = pd.read_csv(os.path.join(root_dir, "metadata/annotation.all.csv"))

# get guide annotation
guide_annotation = os.path.join(root_dir, "metadata/guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)

# # merge output (reads in constructs and assignemnts) of each sample
# reads_all = scores_all = assignment_all = pd.DataFrame()

# for merged_sample_name in sample_annotation["merged_sample_name"].unique():
#     reads = scores = assignment = exp = pd.DataFrame()
#     for sample_name in sample_annotation[sample_annotation["merged_sample_name"] == merged_sample_name]["sample_name"].unique():
#         print(merged_sample_name, sample_name)
#         df1 = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "quantification", "guide_cell_quantification.csv"))
#         df2 = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "quantification", "guide_cell_scores.csv"), index_col=0).reset_index()
#         df3 = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "quantification", "guide_cell_assignment.csv"))

#         df1["sample"] = df2["sample"] = df3["sample"] = sample_name
#         df1["experiment"] = df2["experiment"] = df3["experiment"] = merged_sample_name
#         reads = reads.append(df1)
#         scores = scores.append(df2)
#         assignment = assignment.append(df3)

#     reads.to_csv(os.path.join(results_dir, "guide_cell_quantification.{}.csv".format(merged_sample_name)), index=False)
#     scores.to_csv(os.path.join(results_dir, "guide_cell_scores.{}.csv".format(merged_sample_name)), index=False)
#     assignment.to_csv(os.path.join(results_dir, "guide_cell_assignment.{}.csv".format(merged_sample_name)), index=False)

#     reads_all = reads_all.append(reads)
#     scores_all = scores_all.append(scores)
#     assignment_all = assignment_all.append(assignment)

# reads_all.to_csv(os.path.join(results_dir, "guide_cell_quantification.all.csv"), index=False)
# scores_all.to_csv(os.path.join(results_dir, "guide_cell_scores.all.csv"), index=False)
# assignment_all.to_csv(os.path.join(results_dir, "guide_cell_assignment.all.csv"), index=False)

reads = pd.read_csv(os.path.join(results_dir, "guide_cell_quantification.all.csv"))
scores = pd.read_csv(os.path.join(results_dir, "guide_cell_scores.all.csv"))
assignment = pd.read_csv(os.path.join(results_dir, "guide_cell_assignment.all.csv"))

# merge expression
exp_all = pd.DataFrame()
for n_genes in [500, 100]:
    for i, merged_sample_name in enumerate(sample_annotation["merged_sample_name"].unique()):
        for j, sample_name in enumerate(sample_annotation[sample_annotation["merged_sample_name"] == merged_sample_name]["sample_name"].unique()):
            df = pd.read_csv(os.path.join(root_dir, "results_pipeline", sample_name, "digital_expression.{}genes.tsv".format(n_genes)), sep="\t").set_index("GENE")
            print(merged_sample_name, sample_name, df.shape)
            if j == 0:
                exp = df
            else:
                exp = exp.add(df).fillna(0)
        print("saving")
        exp.to_csv(os.path.join(results_dir, "digital_expression.{}genes.{}.csv".format(n_genes, merged_sample_name)), index=True)
        exp.columns = merged_sample_name + ":" + exp.columns
        if i == 0:
            exp_all = exp
        else:
            exp_all = exp_all.join(exp, how="outer")
        print(merged_sample_name, exp_all.shape)
    print("saving big")
    exp_all.to_csv(os.path.join(results_dir, "digital_expression.{}genes.all.csv".format(n_genes)), index=True)


for n_genes in [500, 100]:
    for i, merged_sample_name in enumerate(sample_annotation["merged_sample_name"].unique()):
        exp = pd.read_csv(os.path.join(results_dir, "digital_expression.{}genes.{}.csv".format(n_genes, merged_sample_name)), index_col=0)
        print(merged_sample_name, exp.shape)

        # stats
        cells_per_gene = exp.apply(lambda x: (x > 0).sum(), axis=1)
        genes_per_cell = exp.apply(lambda x: (x > 0).sum(), axis=0)

        # plot number of genes assigned depending on the minimum number of genes required
        total = list()
        assigned = list()
        for j in range(1, 10000, 10):
            print(j)
            pos = genes_per_cell > j
            d = genes_per_cell[pos]
            try:
                original_names = pd.Series(d.index.str.split("_")).apply(lambda x: x[0])
                assigned.append(original_names.isin(scores["cell"].dropna()).sum())
                total.append(pos.sum())
            except AttributeError:
                break

        fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
        axis[0].plot((range(1, len(total) * 10, 10)), (total), lw=4, color="b", label="All cells")
        axis[0].plot((range(1, len(total) * 10, 10)), (assigned), lw=4, color="r", label="Cells with gRNA assigned")
        axis[1].plot(np.log2(range(1, len(total) * 10, 10)), np.log2(total), lw=4, color="b", label="All cells")
        axis[1].plot(np.log2(range(1, len(total) * 10, 10)), np.log2(assigned), lw=4, color="r", label="Cells with gRNA assigned")
        axis[0].set_xlabel("Genes detected")
        axis[1].set_xlabel("Genes detected (log2)")
        axis[0].set_ylabel("Cells")
        axis[1].set_ylabel("Cells (log2)")
        plt.legend()
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "barcodes_per_cell.{}.identified.varying_genes_detected.svg".format(merged_sample_name)), bbox_inches="tight")

        # cumulative sum
        d = pd.DataFrame([total, assigned], index=['total', 'assigned'], columns=range(1, len(total) * 10, 10)).T
        cumsum = d.cumsum() / d.sum()

        fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
        axis[0].plot(cumsum.index, cumsum['total'], lw=4, color="b", label="All cells")
        axis[0].plot(cumsum.index, cumsum['assigned'], lw=4, color="r", label="Cells with gRNA assigned")
        axis[1].plot(np.log2(cumsum.index), cumsum['total'], lw=4, color="b", label="All cells")
        axis[1].plot(np.log2(cumsum.index), cumsum['assigned'], lw=4, color="r", label="Cells with gRNA assigned")
        axis[0].set_xlabel("Genes detected")
        axis[1].set_xlabel("Genes detected (log2)")
        axis[0].set_ylabel("Cells (fraction of total)")
        axis[1].set_ylabel("Cells (fraction of total)")
        plt.legend()
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "barcodes_per_cell.{}.identified.varying_genes_detected.cumulative.svg".format(merged_sample_name)), bbox_inches="tight")

        # Plot reads vs genes covered for each cell
        reads_per_cell = exp.sum(0)

        fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
        axis[0].scatter(reads_per_cell, genes_per_cell, s=4)
        axis[1].scatter(np.log2(reads_per_cell), np.log2(genes_per_cell), s=4)
        axis[0].set_xlabel("Reads")
        axis[1].set_xlabel("Reads (log2)")
        axis[0].set_ylabel("Genes")
        axis[1].set_ylabel("Genes (log2)")
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "reads_vs_genes.{}.svg".format(merged_sample_name)), bbox_inches="tight")



# exp = pd.read_csv(os.path.join(results_dir, "digital_expression.500genes.HEK293T_4_WNT_unstimulated.csv"), index_col=0)
# exp2 = pd.read_csv(os.path.join(results_dir, "digital_expression.500genes.HEK293T_4_WNT_stimulated.csv"), index_col=0)
# exp = exp.join(exp2)

# plot
# plot_reads_in_constructs(reads)

# Read in expression
# write assignement in name of cell in expression matrix


# Choose thresehold and output matrix
matrix = exp[cells_per_gene > 3][exp.columns[genes_per_cell >= 100]]
matrix.to_csv(os.path.join(root_dir, "digital_expression.3reads_per_gene-3cells_per_gene-100genes_per_cell.tsv"), index=True)

matrix = pd.read_csv(os.path.join(root_dir, "digital_expression.3reads_per_gene-3cells_per_gene-100genes_per_cell.tsv"), index_col=0, sep="\t")

# compare with numbers from summary file
summary = pd.read_csv(os.path.join(root_dir, "digital_expression.summary.100genes.tsv"), sep="\t", skiprows=2, index_col=0)
reads_per_cell = exp.sum(0)
df = pd.DataFrame([reads_per_cell, summary['NUM_GENES']], index=["reads", "genes"]).T

fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
axis[0].scatter(df['reads'], df['genes'], s=4)
axis[1].scatter(np.log2(df['reads']), np.log2(df['genes']), s=4)
axis[0].set_xlabel("Reads")
axis[1].set_xlabel("Reads (log2)")
axis[0].set_ylabel("Genes")
axis[1].set_ylabel("Genes (log2)")
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "reads_vs_genes.from_summary.svg"), bbox_inches="tight")

# select genes present in at least n cells
# select cells with at least n genes
for n_cells in [3, 10, 100]:
    for n_genes in [100, 500, 1000]:
        # n_cells = 100
        # n_genes = 100
        matrix = exp[cells_per_gene > n_cells][exp.columns[genes_per_cell >= n_genes]]

        # transform
        matrix = np.log2(1 + matrix)

        # normalize by total
        matrix_norm = matrix.apply(lambda x: x / x.sum(), axis=0) * 1e4

        # expression heatmaps
        # all genes/cells
        # f = big_heatmap(np.log2(1 + exp).apply(lambda x: x / x.sum(), axis=0) * 1e4, assignment)
        # f.savefig(os.path.join(results_dir, "expression_matrix.total.png"), bbox_inches="tight", dpi=300)
        # reduced
        for name, m in [("", matrix), ("norm-", matrix_norm)]:
            f = big_heatmap(m, assignment)
            f.savefig(os.path.join(results_dir, "expression_matrix.%s%icells_%igenes.png" % (name, n_cells, n_genes)), bbox_inches="tight", dpi=300)

            f = big_heatmap(m[m.columns[m.columns.isin(assignment['cell'])]], assignment)
            f.savefig(os.path.join(results_dir, "expression_matrix.%s%icells_%igenes.only_assigned.png" % (name, n_cells, n_genes)), bbox_inches="tight", dpi=300)

        # Inspect
        from sklearn.decomposition import PCA
        from sklearn.manifold import TSNE, MDS, LocallyLinearEmbedding, SpectralEmbedding, Isomap
        from matplotlib import cm

        # color mapping
        integer_map = dict([(val, i) for i, val in enumerate(set(pd.Series(exp.columns.str.split("_")).apply(lambda x: x[1] if len(x) > 1 else "nan")))])
        colors = [cm.viridis(integer_map[x[1]] / 4.) if len(x) > 1 else cm.viridis(integer_map["nan"] / 4.) for x in matrix.columns.str.split("_")]

        fig, axis = plt.subplots(3, 2, figsize=(12, 8))
        axis = axis.flatten()
        for i, method in enumerate([PCA, MDS, TSNE, LocallyLinearEmbedding, SpectralEmbedding, Isomap]):
            fitted = method().fit_transform(matrix.T)
            axis[i].scatter(fitted[:, 0], fitted[:, 1], color=colors)
        fig.savefig(os.path.join(results_dir, "clustering.%icells_%igenes.png" % (n_cells, n_genes)), bbox_inches="tight")

        fig, axis = plt.subplots(3, 2, figsize=(12, 8))
        axis = axis.flatten()
        for i, method in enumerate([PCA, MDS, TSNE, LocallyLinearEmbedding, SpectralEmbedding, Isomap]):
            fitted = method().fit_transform(matrix_norm.T)
            axis[i].scatter(fitted[:, 0], fitted[:, 1], color=colors)
        fig.savefig(os.path.join(results_dir, "clustering.norm-%icells_%igenes.png" % (n_cells, n_genes)), bbox_inches="tight")

sns.clustermap(matrix_norm.ix[['TET2_gene', 'MBD1_gene', 'filler_gene']])

sns.clustermap(matrix_norm.ix[['TET2', 'TET2_gene', 'MBD1_gene', 'filler_gene']])


# Plot
reads = pd.read_csv(os.path.join(results_dir, "guide_cell_quantification.csv"))
scores = pd.read_csv(os.path.join(results_dir, "guide_cell_scores.csv"), index_col=0)
assignment = pd.read_csv(os.path.join(results_dir, "guide_cell_assignment.csv"))

# calculate abs amount of basepairs overlaping the gRNA of the assigned cell vs all others
overlap_per_cell = reads.groupby(["cell"])['overlap'].sum()
overlap_per_guide = reads.groupby(["cell", "chrom"])['overlap'].sum()
overlap_assignment = scores.apply(lambda x: overlap_per_guide.ix[x.name, x['assignment']] if not pd.isnull(x['assignment']) else pd.np.nan, axis=1)
overlap_others = overlap_per_cell - overlap_assignment

sns.jointplot(overlap_others, overlap_assignment, alpha=0.1)
plt.savefig(os.path.join(results_dir, "duplets_assignment_overlap.svg"), bbox_inches="tight")
plt.close("all")
sns.jointplot(overlap_others, np.log2(1 + overlap_assignment), alpha=0.1)
plt.savefig(os.path.join(results_dir, "duplets_assignment_overlap.ylog.svg"), bbox_inches="tight")
plt.close("all")
sns.jointplot(np.log2(1 + overlap_others), np.log2(1 + overlap_assignment), alpha=0.1)
plt.savefig(os.path.join(results_dir, "duplets_assignment_overlap.bothlog.svg"), bbox_inches="tight")
plt.close("all")
sns.jointplot(overlap_others, overlap_assignment, xlim=(-100, overlap_assignment.max() + 100), ylim=(-100, overlap_assignment.max() + 100), alpha=0.1)
plt.savefig(os.path.join(results_dir, "duplets_assignment_overlap.lims.svg"), bbox_inches="tight")
plt.close("all")


# Count duplets
overlap_others = overlap_others.dropna()
print((overlap_others[overlap_others > 10].shape[0] / float(len(reads['cell'].unique()))) * 100)

# Scatter plots between gRNAs
matrix = np.log2(1 + exp[cells_per_gene > 3][exp.columns[genes_per_cell >= 100]]).apply(lambda x: x / x.sum(), axis=0) * 1e4

# Get "differential"
diff = np.log2(
    (1 + matrix[[x for x in assignment[assignment['assignment'] == "MBD1"]['cell'].tolist() if x in matrix.columns]].mean(1)) /
    (1 + matrix[[x for x in assignment[assignment['assignment'] == "filler"]['cell'].tolist() if x in matrix.columns]].mean(1))).sort_values()

fig, axis = plt.subplots(2, 2, figsize=(12, 8))
axis = iter(axis.flatten())
for i, gene in enumerate(["MBD1", "TET2", "DNMT3B"]):
    i = (1 + matrix[[x for x in assignment[assignment['assignment'] == gene]['cell'].tolist() if x in matrix.columns]].mean(1))
    j = (1 + matrix[[x for x in assignment[assignment['assignment'] == "filler"]['cell'].tolist() if x in matrix.columns]].mean(1))
    a = axis.next()
    a.scatter(i, j)
    a.set_title(gene)
for a in axis:
    a.set_axis_off()
sns.despine(fig)
fig.savefig(os.path.join(results_dir, "guide_scatter.mean.png"), bbox_inches="tight", dpi=300)

fig, axis = plt.subplots(1, figsize=(12, 8))
sns.distplot(matrix[[x for x in assignment[assignment['assignment'] == "MBD1"]['cell'].tolist() if x in matrix.columns]].ix['HIST1H4C'], label="MBD1 gRNA", ax=axis, bins=100)
sns.distplot(matrix[[x for x in assignment[assignment['assignment'] == "filler"]['cell'].tolist() if x in matrix.columns]].ix['HIST1H4C'], label="filler gRNA", ax=axis, bins=100)
axis.legend()
fig.savefig(os.path.join(results_dir, "guide_scatter.HIST1H4C_distibution.png"), bbox_inches="tight", dpi=300)


# Investigate duplicates
df = pd.read_csv("cell_umi_barcodes.100genes.tsv", sep="\t")
cd19neg = pd.read_csv("../DROPseq_CLL4_pool1_CD19neg_16cycles/cell_umi_barcodes.100genes.tsv", sep="\t")
cd19pos13 = pd.read_csv("../DROPseq_CLL4_pool1_CD19pos_13cycles/cell_umi_barcodes.100genes.tsv", sep="\t")
cd19pos16 = pd.read_csv("../DROPseq_CLL4_pool1_CD19pos_16cycles/cell_umi_barcodes.100genes.tsv", sep="\t")

fig, axis = plt.subplots(2, 2, figsize=(12, 12), sharey=False)
sns.distplot(np.log2(1 + df["Num_Obs"]), kde=False, bins=100, ax=axis[0][0])
sns.distplot(np.log2(1 + cd19neg["Num_Obs"]), kde=False, bins=100, ax=axis[0][1])
sns.distplot(np.log2(1 + cd19pos13["Num_Obs"]), kde=False, bins=100, ax=axis[1][0])
sns.distplot(np.log2(1 + cd19pos16["Num_Obs"]), kde=False, bins=100, ax=axis[1][1])
axis[0][0].set_title("CROP-seq resequencing - unique rate {}".format((df["Num_Obs"] == 1).sum() / float((df["Num_Obs"] != 1).sum())))
axis[0][1].set_title("CLL CD19 neg 16 cycles - unique rate {}".format((cd19neg["Num_Obs"] == 1).sum() / float((cd19neg["Num_Obs"] != 1).sum())))
axis[1][0].set_title("CLL CD19 pos 13 cycles - unique rate {}".format((cd19pos13["Num_Obs"] == 1).sum() / float((cd19pos13["Num_Obs"] != 1).sum())))
axis[1][1].set_title("CLL CD19 pos 16 cycles - unique rate {}".format((cd19pos16["Num_Obs"] == 1).sum() / float((cd19pos16["Num_Obs"] != 1).sum())))
fig.savefig(os.path.join(results_dir, "duplicate_read.comparison.svg"), bbox_inches="tight")


# Heatmap of WNt genes
matrix = np.log2(1 + matrix)
matrix_norm = matrix.apply(lambda x: x / x.sum(), axis=0) * 1e4

wnt = pd.read_csv("metadata/wnt_qPCR_assay_genes.csv", header=None).squeeze()
sns.clustermap(matrix_norm.ix[wnt].dropna().T, yticklabels=False, figsize=(18, 12))
plt.savefig(os.path.join(results_dir, "wnt_pathway.heatmap.png"), bbox_inches="tight", dpi=300)
