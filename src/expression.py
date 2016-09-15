#!/usr/bin/env python

import argparse
import os
import pandas as pd
import pysam
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


def normalize(df, experiment="", kind="total"):
    def normalize_by_total(df):
        """
        Normalize expression by total number of counts per cell.
        """
        return df.apply(lambda x: (x / x.sum()) * 1e4, axis=0)

    def seurat_regress(df):
        """
        """
        import rpy2.robjects as robj
        import pandas.rpy.common as com

        run = robj.r("""
            function(to_norm, output){
                library(Seurat)
                library(dplyr)
                library(Matrix)

                to_norm = Matrix(as.matrix(to_norm), sparse = TRUE)

                # Initialize the Seurat object with the raw data
                exp <- new("seurat", raw.data = to_norm)

                # Keep all genes expressed in 1% of cells, keep all cells with >= 200 genes
                exp <- Setup(
                    exp, min.cells=dim(to_norm)[2] / 100., min.genes=200,
                    do.logNormalize=T, total.expr=1e4, project="cropseq")

                # Calculate the percentage of mitochondrial genes and store it.
                mito.genes <- grep("^MT-", rownames(exp@data), value=T)
                percent.mito <- colSums(expm1(exp@data[mito.genes, ]))/colSums(expm1(exp@data))

                # AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
                exp <- AddMetaData(exp, percent.mito, "percent.mito")

                # Filter out cells with more than the 95th percentil of mitochondrial transcripts
                exp_subset <- SubsetData(
                    exp, subset.name="percent.mito",
                    accept.high=quantile(percent.mito, 0.95))

                # Regress out
                exp_subset_regressed <- RegressOut(exp_subset, latent.vars=c("nUMI", "percent.mito"))

                write.table(as.matrix(exp_subset_regressed@data), file=output, sep=",", quote=FALSE)

                return(as.matrix(exp_subset_regressed@data))
            }

        """)

        # convert to Python objects
        norm = com.convert_robj(
            run(df,
                os.path.join(results_dir, "digital_expression.500genes.{}.seurat_regressed.csv".format(experiment))))

        return norm

    if kind == "total":
        norm = normalize_by_total(df)
        norm.to_csv(os.path.join(results_dir, "digital_expression.500genes.{}.seurat_regressed.csv".format(experiment)))
        return norm
    else:
        regressed = seurat_regress(df).to_sparse(fill_value=0)
        regressed.to_pickle(os.path.join(results_dir, "digital_expression.500genes.{}.seurat_regressed.pickle".format(experiment)))
        regressed.to_csv(os.path.join(results_dir, "digital_expression.500genes.{}.seurat_regressed.csv".format(experiment)))
        return regressed


def unsupervised(df, prefix=""):
    # Inspect
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE, MDS, LocallyLinearEmbedding, SpectralEmbedding, Isomap
    from matplotlib import cm

    # color mapping
    integer_map = dict([(val, i) for i, val in enumerate(set(pd.Series(df.columns.str.split("_")).apply(lambda x: x[1] if len(x) > 1 else "nan")))])
    colors = [cm.viridis(integer_map[x[1]] / 4.) if len(x) > 1 else cm.viridis(integer_map["nan"] / 4.) for x in matrix.columns.str.split("_")]

    fig, axis = plt.subplots(3, 2, figsize=(12, 8))
    axis = axis.flatten()
    for i, method in enumerate([PCA, MDS, TSNE, LocallyLinearEmbedding, SpectralEmbedding, Isomap]):
        fitted = method().fit_transform(matrix.T)
        axis[i].scatter(fitted[:, 0], fitted[:, 1], color=colors)
    fig.savefig(os.path.join(results_dir, "clustering.{}.png".format(prefix)), bbox_inches="tight")


def differential_genes(pos, neg, df, assignment, prefix="", method="mannwhitney"):
    """
    """
    from scipy.stats import mannwhitneyu
    import multiprocessing
    import parmap
    from statsmodels.sandbox.stats.multicomp import multipletests

    def mannwhitneyu_test(g):
        """
        Test two groups of cells for differences in gene `g`
        using Mann-Whitney's U test.
        Report also log fold-change.
        """
        a = pos.ix[g]
        b = neg.ix[g]
        a_ = a.mean()
        b_ = b.mean()
        return [a_, b_, np.log2((a_) / (b_))] + list(mannwhitneyu(a, b))

    def test_2():
        return

    print("Doing differnetial gene expression for experiment: '{}'".format(experiment))
    print("Comparing expression of {} with {} cells in {} genes.".format(pos.shape[1], neg.shape[1], df.shape[0]))

    pos["intercept"] = 0.1
    neg["intercept"] = 0.1
    # apply test function (manwhitney-U)
    stats = pd.DataFrame(
        map(
            lambda x:
                pd.Series(x),
                parmap.map(
                    mannwhitneyu_test if method == "mannwhitney" else test_2,
                    pos.index,
                    parallel=True
                )
        )
    )
    stats.columns = ["a_mean", "b_mean", "fold_change", "stat", "p_value"]
    stats.index = pos.index.tolist()
    stats["q_value"] = multipletests(stats["p_value"], method="fdr_bh")[1]
    stats["log_difference"] = stats["a_mean"] - stats["b_mean"]
    stats["log_p_value"] = -np.log10(stats["p_value"])
    stats["log_q_value"] = -np.log10(stats["q_value"])
    stats.index.name = "gene"
    stats.to_csv(os.path.join(results_dir, "differential_expression.{}.stimutation.csv".format(prefix)), index=True)

    return stats


def plot_deg_stats(stats, prefix=""):
    """
    """

    # Volcano plot
    colors = [sns.color_palette("colorblind")[2] if x < 0.05 else sns.color_palette("colorblind")[0] for x in stats["q_value"]]
    sizes = [15 if x < 0.05 else 7.5 for x in stats["q_value"]]

    fig, axis = plt.subplots(1, 3, figsize=(18, 6))
    axis[0].scatter(stats["fold_change"], stats["log_p_value"], color=colors, alpha=0.2, s=sizes)
    axis[0].set_title("Volcano plot for {}".format(experiment))
    axis[0].set_xlabel("Log fold change")
    axis[0].set_xlabel("p value")

    axis[1].scatter(stats["fold_change"], stats["log_q_value"], color=colors, alpha=0.2, s=sizes)
    axis[1].set_title("Volcano plot for {}".format(experiment))
    axis[1].set_xlabel("Log fold change")
    axis[1].set_xlabel("q value (BH)")

    # MA plot
    axis[2].scatter(np.log2(pos.mean(1) * neg.mean(1)) / 2., stats["fold_change"], color=colors, alpha=0.2, s=sizes)
    axis[2].set_title("MA plot for {}".format(experiment))
    axis[2].set_xlabel("M")
    axis[2].set_xlabel("A")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.stimutation.png".format(prefix)), bbox_inches="tight", dpi=300)


def plot_deg_heatmaps(df, assignment, stats, prefix=""):
    """
    """
    from matplotlib.colors import colorConverter
    from sklearn.manifold import MDS, LocallyLinearEmbedding, Isomap, SpectralEmbedding, TSNE
    from sklearn.decomposition import PCA

    def get_grna_colors(d, assignment):
        # make color dict
        pallete = sns.color_palette("colorblind") * 10
        gRNA_color_dict = dict(zip(assignment["group"].unique(), pallete))
        gRNA_color_dict["CTRL"] = "#228B22"
        gRNA_color_dict["Essential"] = "#B22222"
        gRNA_color_dict["Unassigned"] = "grey"

        # match colors with matrix
        stimulation_colors = ["#B22222" if "un" in x else "#228B22" for x in d.columns]
        cell_names = d.columns.str.lstrip("un|st")
        ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
        ass = [x if type(x) == str else "Unassigned" for x in ass]
        gRNA_colors = [gRNA_color_dict[x] for x in ass]

        # convert to rgb
        stimulation_colors = [colorConverter.to_rgb(x) if type(x) is str else x for x in stimulation_colors]
        gRNA_colors = [colorConverter.to_rgb(x) if type(x) is str else x for x in gRNA_colors]

        return [stimulation_colors, gRNA_colors]

    def get_group_colors(d, assignment):
        # make color dict
        pallete = sns.color_palette("colorblind") * 10
        gRNA_color_dict = dict(zip(d.columns.levels[1], pallete))
        gRNA_color_dict["CTRL"] = "#228B22"
        gRNA_color_dict["Essential"] = "#B22222"
        gRNA_color_dict["Unassigned"] = "grey"

        # match colors with matrix
        gRNA_colors = [gRNA_color_dict[x] for x in d.columns.get_level_values(1)]
        stimulation_colors = ["#B22222" if "un" in x else "#228B22" for x in d.columns.get_level_values(0)]

        # convert to rgb
        stimulation_colors = [colorConverter.to_rgb(x) if type(x) is str else x for x in stimulation_colors]
        gRNA_colors = [colorConverter.to_rgb(x) if type(x) is str else x for x in gRNA_colors]

        return [stimulation_colors, gRNA_colors]

    def get_foldchange_colors(d, s):
        fold_change_colors = ["#B22222" if x > 0 else "#333399" for x in s.ix[d.index]["fold_change"]]
        return fold_change_colors

    # get top 500 differential expressed (if significant) in top 2000 cells (by coverage)
    df2 = df.ix[
        stats[
            (stats["q_value"] < 0.05) &  # filter by p-value
            (abs(stats["fold_change"]) > np.log2(1.5))  # filter by fold-change
        ].sort_values("q_value").head(500).index.tolist()
    ]
    # save stats of diff genes
    stats.ix[df2.index.tolist()].to_csv(os.path.join(results_dir, "differential_expression.{}.differential_genes.csv".format(prefix)), index=True)

    # get top 500 cells from each condition
    df3 = df2[
        df2[df2.columns[df2.columns.str.contains("st")]].sum().sort_values().tail(1500).index.tolist() +
        df2[df2.columns[df2.columns.str.contains("un")]].sum().sort_values().tail(1500).index.tolist()
    ]

    # remove opposite gRNA library
    if "TCR" in experiment:
        assignment = assignment[~assignment['assignment'].str.contains("Wnt")]
    elif "WNT" in experiment:
        assignment = assignment[~assignment['assignment'].str.contains("Tcr")]

    g = sns.clustermap(
        df3, z_score=0,
        col_colors=get_grna_colors(df3, assignment),
        xticklabels=False, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.clustermap.png".format(prefix)), bbox_inches="tight", dpi=300)

    # sort rows by fold change, columns by stimulation
    df4 = df3.ix[stats['fold_change'].sort_values().index].dropna()
    df4 = df4[df4.columns[df4.sum() != 0]]
    df4 = df4.T.sort_values(df4.index.tolist(), ascending=False).T

    g = sns.clustermap(
        df4, z_score=1,
        col_colors=get_grna_colors(df4, assignment),
        row_cluster=False, col_cluster=False,
        xticklabels=False, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sorted_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)

    # sort by stimulation order
    df4 = df3.sort_index(axis=1)

    g = sns.clustermap(
        df4, z_score=0,
        col_colors=get_grna_colors(df4, assignment),
        row_cluster=True, col_cluster=False,
        xticklabels=False, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sortedcondition_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)

    #

    # sort by stimulation order and gRNA order
    cell_names = df3.columns.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    ass = [x if type(x) == str else "Unassigned" for x in ass]
    df3.loc["ass", :] = ass
    df3.loc["sti", :] = [x[:2] for x in df3.columns]
    order = df3.T[["sti", "ass"]].sort_values(["sti", "ass"]).index.tolist()
    df3 = df3.drop(["sti", "ass"])
    df3 = df3.astype(np.float64)

    df4 = df3[order]

    g = sns.clustermap(
        df4, z_score=0,
        col_colors=get_grna_colors(df4, assignment),
        row_cluster=True, col_cluster=False,
        xticklabels=False, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sortedconditiongRNA_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)

    #

    # Group by stimulation / gRNA, get mean expression

    # use all cells for this
    df4 = df2.copy().T

    cell_names = df4.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    ass = [x if type(x) == str else "Unassigned" for x in ass]

    df4["ass"] = ass
    df4['ass'] = df4['ass'].astype("category")
    df4["sti"] = [x[:2] for x in df4.index]
    df4['sti'] = df4['sti'].astype("category")

    df5 = df4.groupby(['sti', 'ass']).mean().T

    # report number of cells per condition & gRNA
    fig, axis = plt.subplots(1, figsize=(6, 12))
    c = df4.groupby(['sti', 'ass']).apply(len).sort_values(ascending=False)
    sns.barplot(c, c.index.to_series().str.join(" "), orient="horiz", ax=axis)
    axis.set_xscale('log')
    axis.set_xlabel('Number of cells (log)')
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.cell_number_per_group.svg".format(prefix)), bbox_inches="tight")

    # Filter out genes with less than 10 cells and save
    df5 = df5[c[c >= 10].index]
    df5.to_csv(os.path.join(results_dir, "differential_expression.{}.group_means.csv".format(prefix)), index=True)

    # cluster
    g = sns.clustermap(
        df5, z_score=0,
        col_colors=get_group_colors(df5, assignment),
        row_colors=get_foldchange_colors(df5, stats),
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.clustered.png".format(prefix)), bbox_inches="tight", dpi=300)

    # sort by stimulation order and gRNA order
    g = sns.clustermap(
        df5.T.sort_index().T,  # z_score=0,
        col_colors=get_group_colors(df5.T.sort_index().T, assignment),
        row_colors=get_foldchange_colors(df5, stats),
        row_cluster=True, col_cluster=False,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.sortedconditiongRNA_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)

    #

    # Dimentionality reduction methods
    # try several
    methods = ["PCA", "LocallyLinearEmbedding", "Isomap", "SpectralEmbedding", "TSNE", "MDS"]

    for name, matrix in [("groups", df5.T), ("cells", df3.T)]:
        for method in methods:
            print(name, method)
            model = eval(method)()

            if name == "cells":
                color = get_grna_colors(matrix.T, assignment)
                m = matrix.T[matrix.dtypes == np.float64].T
            else:
                color = get_group_colors(matrix.T, assignment)
                m = matrix.T[matrix.dtypes == np.float64]

            fit = model.fit_transform(m)

            # plot
            if method == "PCA":
                pcs = 3
            else:
                pcs = 1

            fig, axis = plt.subplots(2, pcs, sharex=True, sharey=True, figsize=(8, 10))
            if method != "PCA":
                axis = [[x] for x in axis]

            for i, variable in enumerate(["condition", "gene"]):
                for pc in range(pcs):
                    axis[i][pc].scatter(fit[:, pc], fit[:, pc + 1], color=color[i], alpha=0.75 if name == "groups" else 0.1)
                    axis[i][pc].set_xticklabels([])
                    axis[i][pc].set_yticklabels([])
                    if method == "PCA":
                        axis[i][pc].set_xlabel("PC %i" % pc)
                        axis[i][pc].set_ylabel("PC %i" % (pc + 1))
                # axis[i][pc].legend(
                #     handles=[mpatches.Patch(color=v, label=k) for k, v in color_mapping.items()],
                #     ncol=2 if feature == "organ" else 1,
                #     loc='center left',
                #     bbox_to_anchor=(1, 0.5))
            fig.savefig(os.path.join(results_dir, "differential_expression.{}.{}.{}.png".format(prefix, name, method)), bbox_inches="tight", dpi=300)


def enrich_signature(degs, prefix=""):
    """
    """
    # load stats of diff genes
    degs = pd.read_csv(os.path.join(results_dir, "differential_expression.{}.differential_genes.csv".format(prefix)), index_col=0)
    degs.index.name = "gene_name"

    for d, name in [(degs, "_all_genes"), (degs[degs["fold_change"] > 0], "_up_genes"), (degs[degs["fold_change"] < 0], "_down_genes")]:

        enr = enrichr(d.reset_index())
        enr.to_csv(os.path.join(results_dir, "differential_expression.{}.enrichr.csv".format(prefix + name)), index=False)

        for gene_set_library in enr["gene_set_library"].unique():

            p = enr[enr["gene_set_library"] == gene_set_library]
            # p = p[p["adjusted_p_value"] < 0.05]

            if p.shape[0] < 1:
                continue

            p = p.sort_values("combined_score").tail(25)

            print(gene_set_library)
            fig, axis = plt.subplots(1)
            sns.barplot(p["combined_score"], p["description"], orient="horiz", ax=axis)
            axis.set_xlabel("Combined score")
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "differential_expression.{}.enrichr.{}.svg".format(prefix + name, gene_set_library)), bbox_inches="tight")


def assign_cells_to_signature(degs, df, assignment, prefix=""):
    """
    """
    group_means = pd.read_csv(
        os.path.join(results_dir, "differential_expression.{}.group_means.csv".format(prefix)),
        index_col=0,
        header=[0, 1], skipinitialspace=True, tupleize_cols=True)
    group_means.columns = pd.MultiIndex.from_tuples(group_means.columns)

    # Get trait-specific signature
    # 1. get median accessibility of each group
    x1 = group_means["st", "CTRL"]
    x2 = group_means["un", "CTRL"]

    # 2. get signature matrix
    # here, bounds are set to (-20, 20) so that the extremes of the signature represent -20% or 120% of the signature
    # this is done because the extreme values (0 and 100%) values correspond to the median value within each group,
    # meaning that some samples are expected to surpass those values.
    sign = generate_signature_matrix(np.vstack([x1, x2]).T, n=101, bounds=(-20, 20))

    # 3. get signature value of each patient
    # x[[s.name for s in samples]].apply(best_signature_matrix, matrix=sign, axis=1)
    sigs = list()
    for i, cell in enumerate(df.columns):
        if i % 20 == 0:
            print(i)
        sigs.append(best_signature_matrix(array=df.ix[group_means.index][cell], matrix=sign))

    # Save clinical trait signatures for all samples
    sigs.to_csv(os.path.join(results_dir, "signatures.all_cells.{}.csv".format(prefix)))

    # Aggregate by condition/gene
    df2 = df.copy().T
    # use all cells for this
    cell_names = df2.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    ass = [x if type(x) == str else "Unassigned" for x in ass]

    df2["ass"] = ass
    df2['ass'] = df2['ass'].astype("category")
    df2["sti"] = [x[:2] for x in df2.index]
    df2['sti'] = df2['sti'].astype("category")

    # add signatures
    df2['signature'] = pd.Series(sigs).astype(np.int64)

    # plot distibution
    sigs_mean = df2.groupby(['sti', 'ass'])['signatures'].mean().sort_values()

    fig, axis = plt.subplots(1, figsize=(10, 8))
    sns.stripplot(x=sigs_mean, y=sigs_mean.index, orient="horiz", ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.strip.svg".format(prefix)), bbox_inches="tight")

    # Given the known condition of each group, what is the deviation from that?
    # plot as rank of mean
    fig, axis = plt.subplots(1, figsize=(10, 8))
    axis.scatter(sigs_mean.rank(ascending=False), sigs_mean)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.rank.svg".format(prefix)), bbox_inches="tight")

    # plot as violinplots (values per cell)
    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(x="sti", y="signatures", hue="ass", data=df2, ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.violinplot.svg".format(prefix)), bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(x="ass", y="signatures", hue="sti", data=df2, ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.violinplot2.svg".format(prefix)), bbox_inches="tight")

    # Calculate deviation from CTRL for each stimulation
    s = (sigs_mean.ix['st'] - sigs_mean.ix['st', "CTRL"]).sort_values(ascending=False)
    u = (sigs_mean.ix['un'] - sigs_mean.ix['un', "CTRL"]).sort_values(ascending=False)
    sl = np.log2(sigs_mean.ix['st'] / sigs_mean.ix['st', "CTRL"]).sort_values(ascending=False)
    ul = np.log2(sigs_mean.ix['un'] / sigs_mean.ix['un', "CTRL"]).sort_values(ascending=False)

    fig, axis = plt.subplots(2, figsize=(8, 8))
    sns.barplot(s, s.index, orient="horiz", order=s.index, ax=axis[0])
    sns.barplot(u, u.index, orient="horiz", order=u.index, ax=axis[1])
    axis[0].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_ylabel("Stimulated")
    axis[1].set_ylabel("Unstimulated")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature_deviation.rank.barplot.svg".format(prefix)), bbox_inches="tight")

    fig, axis = plt.subplots(2, figsize=(8, 8))
    sns.barplot(sl, sl.index, orient="horiz", order=sl.index, ax=axis[0])
    sns.barplot(ul, ul.index, orient="horiz", order=ul.index, ax=axis[1])
    axis[0].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_ylabel("Stimulated")
    axis[1].set_ylabel("Unstimulated")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature_deviation.ranklog2.barplot.svg".format(prefix)), bbox_inches="tight")

    fig, axis = plt.subplots(2, figsize=(8, 8))
    axis[0].scatter(s.rank(ascending=False), s)
    axis[1].scatter(u.rank(ascending=False), u)
    axis[0].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_title("Stimulated")
    axis[1].set_title("Unstimulated")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature_deviation.rank.scatter.svg".format(prefix)), bbox_inches="tight")

    fig, axis = plt.subplots(2, figsize=(8, 8))
    axis[0].scatter(sl.rank(ascending=False), sl)
    axis[1].scatter(ul.rank(ascending=False), ul)
    axis[0].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_title("Stimulated")
    axis[1].set_title("Unstimulated")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature_deviation.ranklog2.scatter.svg".format(prefix)), bbox_inches="tight")


def generate_signature_matrix(array, n=101, bounds=(0, 0)):
    """
    :param np.array: 2D np.array
    """
    def get_score(i, j, p, n):
        """Get signature score between p% of the values of the two groups."""
        return ((float(i) * p) + (float(j) * (n - p))) / n

    matrix = np.zeros([array.shape[0], n])
    for x in range(array.shape[0]):
        for y, p in enumerate(np.linspace(0 + bounds[0], n + bounds[1], n)):
            matrix[x, y] = get_score(array[x, 0], array[x, 1], p, n)

    return matrix


def best_signature_matrix(array, matrix):
    """
    :param np.array: 2D np.array
    """
    from scipy.stats import pearsonr

    cors = dict()
    for i in range(matrix.shape[1]):
        cors[i] = pearsonr(array, matrix[:, i])

    return cors.values().index(max(cors.values()))  # index
    # (
    # cors.values().index(max(cors.values())),  # index
    # max(cors.values())  # highest correlation value


def enrichr(dataframe, gene_set_libraries=None, kind="genes"):
    """
    Use Enrichr on a list of genes (currently only genes supported through the API).
    """
    import json
    import requests

    ENRICHR_ADD = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    ENRICHR_RETRIEVE = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    if gene_set_libraries is None:
        gene_set_libraries = [
            'GO_Biological_Process_2015',
            'GO_Molecular_Function_2015',
            'GO_Cellular_Component_2015',
            "ChEA_2015",
            "KEGG_2016",
            "WikiPathways_2016",
            "Reactome_2016",
            "BioCarta_2016",
            "Disease_Perturbations_from_GEO_down",
            "Disease_Perturbations_from_GEO_up",
            "TF-LOF_Expression_from_GEO"
        ]

    results = pd.DataFrame()
    for gene_set_library in gene_set_libraries:
        print("Using enricher on %s gene set library." % gene_set_library)

        if kind == "genes":
            # Build payload with bed file
            attr = "\n".join(dataframe["gene_name"].dropna().tolist())
        elif kind == "regions":
            # Build payload with bed file
            attr = "\n".join(dataframe[['chrom', 'start', 'end']].apply(lambda x: "\t".join([str(i) for i in x]), axis=1).tolist())

        payload = {
            'list': (None, attr),
            'description': (None, gene_set_library)
        }
        # Request adding gene set
        response = requests.post(ENRICHR_ADD, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        # Track gene set ID
        user_list_id = json.loads(response.text)['userListId']

        # Request enriched sets in gene set
        response = requests.get(
            ENRICHR_RETRIEVE + query_string % (user_list_id, gene_set_library)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        # Get enriched sets in gene set
        res = json.loads(response.text)
        # If there's no enrichemnt, continue
        if len(res) < 0:
            continue

        # Put in dataframe
        res = pd.DataFrame([pd.Series(s) for s in res[gene_set_library]])
        if res.shape[0] == 0:
            continue
        res.columns = ["rank", "description", "p_value", "z_score", "combined_score", "genes", "adjusted_p_value"]

        # Remember gene set library used
        res["gene_set_library"] = gene_set_library

        # Append to master dataframe
        results = results.append(res, ignore_index=True)

    return results


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
sample_annotation = pd.read_csv(os.path.join(root_dir, "metadata/annotation.csv"))

# get guide annotation
guide_annotation = os.path.join(root_dir, "metadata/guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)

reads = pd.read_csv(os.path.join(results_dir, "guide_cell_quantification.all.csv"))
scores = pd.read_csv(os.path.join(results_dir, "guide_cell_scores.all.csv"))
assignment = pd.read_csv(os.path.join(results_dir, "guide_cell_assignment.all.csv"))

# group gRNAs per gene or control group
assignment = assignment.set_index("assignment")
assignment["group"] = pd.np.nan
assignment.loc[assignment.index.str.contains("CTRL"), "group"] = "CTRL"
assignment.loc[assignment.index.str.contains("Essential"), "group"] = "Essential"
assignment = assignment.reset_index()
tr = assignment[assignment["group"].isnull()]['assignment'].str.extract(".*_(.*)_.*")
assignment.loc[tr.index, "group"] = tr

# group by experiment (without stimulation condition)
assignment["experiment_group"] = assignment["experiment"].str.extract("^(.*)_.*$")

n_genes = 500
experiment = "CROP-seq_Jurkat_TCR"
prefix = experiment + "_stimulation.allcells"
method = "mannwhitney"

# get expression
for n_genes in [500]:
    for experiment in ["CROP-seq_Jurkat_TCR", "CROP-seq_HEK293T_4_WNT", "CROP-seq_HEK293T_6_WNT"]:
        print(experiment)

        # exp = pd.read_csv(os.path.join(results_dir, "digital_expression.{}genes.{}.csv".format(n_genes, experiment)), index=True).to_sparse(fill_value=0)
        exp = pd.read_pickle(os.path.join(results_dir, "digital_expression.{}genes.{}.pickle".format(n_genes, experiment)))

        # Filter matrix
        cells_per_gene = exp.apply(lambda x: (x > 2).sum(), axis=1)
        genes_per_cell = exp.apply(lambda x: (x > 2).sum(), axis=0)
        matrix = exp[cells_per_gene > 20]
        matrix = matrix[exp.columns[genes_per_cell >= 100]]

        #

        # Normalize

        # Approach 1:
        # normalize by total
        matrix_norm = normalize(matrix, experiment=experiment, kind="total")

        # Approach 2:
        # regress out based on total number and MT-genes using Seurat
        matrix_norm = normalize(matrix, kind="seurat")
        matrix_norm = pd.read_csv(os.path.join(results_dir, "digital_expression.500genes.{}.seurat_regressed.csv".format(experiment)), index_col=0).to_sparse(fill_value=0)
        matrix_norm.to_pickle(os.path.join(results_dir, "digital_expression.500genes.{}.seurat_regressed.pickle".format(experiment)))
        matrix_norm = pd.read_pickle(os.path.join(results_dir, "digital_expression.500genes.{}.seurat_regressed.pickle".format(experiment)))

        #

        #

        # Unsupervised

        # Approach 1:
        # apply dimentionality reduction methods/clustering
        # and discover biological axis related with stimulation

        # Approach 2 (just a bit supervised though):
        # get differential genes between conditions from CTRL cells only
        # use this signature to position each cell
        # observe deviation of groups of cells with gRNA targeting the same
        # gene coming from either condition to the mean

        # variant A:
        # all cells from each group
        pos = matrix_norm[matrix_norm.columns[matrix_norm.columns.str.contains("st")].tolist()]
        neg = matrix_norm[matrix_norm.columns[matrix_norm.columns.str.contains("un")].tolist()]

        prefix = experiment + "_stimulation.allcells"

        s = os.path.join(results_dir, "differential_expression.{}.stimutation.csv".format(prefix))
        if not os.path.exists(s):
            stats = differential_genes(pos, neg, matrix_norm, assignment, prefix=prefix)
        else:
            stats = pd.read_csv(s, index_col=0)

        if not os.path.exists(os.path.join(results_dir, "differential_expression.{}.stimutation.png".format(prefix))):
            plot_deg_stats(stats, prefix=prefix)

        if not os.path.exists(os.path.join(results_dir, "differential_expression.{}.group_means.png".format(prefix))):
            plot_deg_heatmaps(matrix_norm.to_dense(), assignment, stats, prefix=prefix)

        if not os.path.join(results_dir, "differential_expression.{}.enrichr.csv".format(prefix + "_all_genes")):
            enrich_signature(stats, prefix=prefix)




        # variant B:
        # only "CTRL" cells
        experiment_assignment = assignment[assignment["experiment_group"] == experiment]
        ctrl_cells = experiment_assignment[experiment_assignment['group'] == "CTRL"]['cell']
        cell_names = matrix_norm.columns.str.lstrip("un|st")

        selected_cells = pd.Series(matrix_norm.columns).where(cell_names.isin(ctrl_cells)).dropna().tolist()

        ctrl_matrix = matrix_norm[selected_cells]

        pos = ctrl_matrix[ctrl_matrix.columns[ctrl_matrix.columns.str.contains("st")].tolist()].to_dense()
        neg = ctrl_matrix[ctrl_matrix.columns[ctrl_matrix.columns.str.contains("un")].tolist()].to_dense()

        stats = differential_genes(pos, neg, ctrl_matrix, assignment, prefix=experiment + "_stimulation.ctrlcells")

        # Supervised
        # Approach 1:
        # Compare groups of cells with gRNA targeting the same gene with either:
        # a) control cells; b) all cells; - intersect (or max) the differential genes from either,
        # explore signature enrichment

        # Approach 2:
        # get MSigDB/GEO signatures on these stimuli
        # use them to quantify the deviation of each group of cell with gRNA targeting the same gene
