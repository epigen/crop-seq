#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from looper.models import Project


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
            module unload R
            module load gcc/6.0.0
            module load texlive
            module load R/3.2.3

            function(to_norm, output){
                library(Seurat)
                library(dplyr)
                library(Matrix)

                # if run manually:
                library(rhdf5)
                to_norm <- h5read("~/projects/crop-seq/results/CROP-seq_Jurkat_TCR.digital_expression.500genes.hdf5", "exp_matrix")

                # else:
                # to_norm = Matrix(as.matrix(to_norm), sparse = TRUE)

                # Make dataframe
                M = as.data.frame(t(to_norm$block0_values))
                names = paste(  # VERY IMPORTANT!: remember R is 1-based and the level indexes come from Python!
                    as.character(to_norm$block0_items_level0[to_norm$block0_items_label0 + 1]),
                    as.character(to_norm$block0_items_level1[to_norm$block0_items_label1 + 1]),
                    as.character(to_norm$block0_items_level2[to_norm$block0_items_label2 + 1])
                    as.character(to_norm$block0_items_level3[to_norm$block0_items_label3 + 1]), sep="")
                colnames(M) <- names
                rownames(M) <- to_norm$axis1

                # Initialize the Seurat object with the raw data
                exp <- new("seurat", raw.data = M)

                # Keep all genes expressed in at least 20 cells, keep all cells with >= 200 genes
                exp <- Setup(
                    exp, project="cropseq",
                    min.cells=20, min.genes=200,
                    do.logNormalize=T, total.expr=1e4)

                # Stash cell attributes
                exp <- AddMetaData(exp, to_norm$block0_items_level0[to_norm$block0_items_label0 + 1], "condition")
                exp <- AddMetaData(exp, as.factor(to_norm$block0_items_level1[to_norm$block0_items_label1 + 1]), "replicate")
                exp <- AddMetaData(exp, to_norm$block0_items_level2[to_norm$block0_items_label2 + 1], "cell")

                # Calculate the percentage of mitochondrial genes and store it.
                mito.genes <- grep("^MT-", rownames(exp@data), value=T)
                percent.mito <- colSums(expm1(exp@data[mito.genes, ])) / colSums(expm1(exp@data))
                exp <- AddMetaData(exp, percent.mito, "percent.mito")
                # Calculate the percentage of ribosomal genes and store it.
                ribo.genes <- grep("^RP", rownames(exp@data), value=T)
                percent.ribo <- colSums(expm1(exp@data[ribo.genes, ])) / colSums(expm1(exp@data))
                exp <- AddMetaData(exp, percent.ribo, "percent.ribo")

                # Plot QC stats
                VlnPlot(exp, c("replicate", "nGene", "nUMI", "percent.mito"))
                GenePlot(exp,"IL2RA","MALAT1",cex.use = 1)
                par(mfrow = c(2, 2))
                GenePlot(exp, "nUMI", "nGene")
                GenePlot(exp, "nUMI", "percent.mito")
                GenePlot(exp, "nUMI", "percent.ribo")
                GenePlot(exp, "percent.ribo", "percent.mito")

                # Filter out cells with more than the 95th percentile of mitochondrial transcripts
                exp_subset <- SubsetData(
                    exp, subset.name="percent.mito",
                    accept.high=quantile(percent.mito, 0.95))

                # Regress out
                exp_subset_regressed <- RegressOut(exp_subset, latent.vars=c("nUMI", "percent.mito"))

                # Write as hdf5
                output = "~/projects/crop-seq/results/CROP-seq_Jurkat_TCR.digital_expression.500genes.seurat_regressed.h5"
                h5createFile(output)
                h5createGroup(output, "seurat_matrix")
                h5write(as.matrix(exp_subset_regressed@scale.data), file=output, "seurat_matrix/matrix")
                h5write(colnames(exp_subset_regressed@scale.data), file=output, "seurat_matrix/columns")
                h5write(rownames(exp_subset_regressed@scale.data), file=output, "seurat_matrix/rows")
                H5close(output)

                write.table(as.matrix(exp_subset_regressed@scale.data), file=output, sep=",", quote=FALSE)
                save(exp_subset_regressed, "~/projects/crop-seq/results/CROP-seq_Jurkat_TCR.digital_expression.500genes.seurat.regressed.Rdata")
                return(as.matrix(exp_subset_regressed@data))


                # continue to explore
                exp_subset_regressed <- MeanVarPlot(exp_subset_regressed ,fxn.x=expMean, fxn.y=logVarDivMean, x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5, do.contour=F)

                length(exp_subset_regressed@var.genes)

                exp_subset_regressed <- PCA(exp_subset_regressed, pc.genes = exp_subset_regressed@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)


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


def read_seurat_hdf5(hdf5_file, assignment):
    import h5py
    with h5py.File(hdf5_file, 'r') as handle:
        cols = handle.get("seurat_matrix/columns").value
        rows = handle.get("seurat_matrix/rows").value
        df = handle.get("seurat_matrix/matrix").value
    seurat_matrix = pd.DataFrame(df, index=cols, columns=rows).T

    assignment['fullname'] = assignment['condition'].astype(str) + assignment['replicate'].astype(int).astype(str) + assignment['cell'].astype(str)

    # add info as multiindex columns
    arrays = [
        seurat_matrix.columns.str.extract("(^.*)\d.*"),
        seurat_matrix.columns.str.extract(".*(\d).*"),
        seurat_matrix.columns.str.extract(".*\d(.*)"),
        assignment.set_index('fullname').ix[seurat_matrix.columns]['assignment'].tolist()]
    seurat_matrix.columns = pd.MultiIndex.from_tuples(list(zip(*arrays)), names=['condition', 'replicate', 'cell', 'grna'])

    return seurat_matrix


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


def differential_genes(pos, neg, assignment, prefix="", method="mannwhitney"):
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

    def scde():
        return

    def deseq2(df):
        df.to_csv("count_matrix.csv")
        pd.DataFrame([df.columns, df.columns.str.get(0)], index=["cell", "condition"]).T.to_csv("col_matrix.csv", index=False)
        return

    print("Doing differnetial gene expression for experiment: '{}'".format(experiment))
    print("Comparing expression of {} with {} cells in {} genes.".format(pos.shape[1], neg.shape[1], pos.shape[0]))

    if method == "deseq":
        return deseq2()

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
    if "WNT" in prefix:
        df2 = stats[(stats["q_value"] < 0.05)]
    else:
        df2 = stats[
            (stats["q_value"] < 0.05) &  # filter by p-value
            (abs(stats["fold_change"]) > np.log2(1.5))  # filter by fold-change
        ]
    sig_down = df2[df2["fold_change"] < 0]
    sig_up = df2[df2["fold_change"] > 0]

    fig, axis = plt.subplots(2, 2, figsize=(12, 12))
    axis = axis.flatten()
    axis[0].hexbin(stats["b_mean"], stats["a_mean"], bins=None, mincnt=1, color="black", alpha=0.2)
    axis[0].scatter(sig_down["b_mean"], sig_down["a_mean"], color=sns.color_palette("colorblind")[0], alpha=0.5, s=10)
    axis[0].scatter(sig_down["b_mean"], sig_down["a_mean"], color=sns.color_palette("colorblind")[0], alpha=0.5, s=10)
    axis[0].scatter(sig_up["b_mean"], sig_up["a_mean"], color=sns.color_palette("colorblind")[2], alpha=0.5, s=10)
    axis[0].set_title("Scatter plot for {}".format(experiment))
    axis[0].set_xlabel(u"Naïve")
    axis[0].set_ylabel("Stimulated")

    axis[1].scatter(stats["fold_change"], stats["log_p_value"], color="grey", alpha=0.2)
    axis[1].scatter(sig_down["fold_change"], sig_down["log_p_value"], color=sns.color_palette("colorblind")[0], alpha=0.2)
    axis[1].scatter(sig_up["fold_change"], sig_up["log_p_value"], color=sns.color_palette("colorblind")[2], alpha=0.2)
    axis[1].set_title("Volcano plot for {}".format(experiment))
    axis[1].set_xlabel("Log fold change")
    axis[1].set_ylabel("p value")

    axis[2].scatter(stats["fold_change"], stats["log_p_value"], color="grey", alpha=0.2)
    axis[2].scatter(sig_down["fold_change"], sig_down["log_p_value"], color=sns.color_palette("colorblind")[0], alpha=0.2)
    axis[2].scatter(sig_up["fold_change"], sig_up["log_p_value"], color=sns.color_palette("colorblind")[2], alpha=0.2)
    axis[2].set_title("Volcano plot for {}".format(experiment))
    axis[2].set_xlabel("Log fold change")
    axis[2].set_ylabel("p value")

    # MA plot
    axis[3].scatter(np.log2(stats["b_mean"] * stats["a_mean"]) / 2., stats["fold_change"], color="grey", alpha=0.2)
    axis[3].scatter(np.log2(sig_down["b_mean"] * sig_down["a_mean"]) / 2., sig_down["fold_change"], color=sns.color_palette("colorblind")[0], alpha=0.2)
    axis[3].scatter(np.log2(sig_up["b_mean"] * sig_up["a_mean"]) / 2., sig_up["fold_change"], color=sns.color_palette("colorblind")[2], alpha=0.2)
    axis[3].set_title("MA plot for {}".format(experiment))
    axis[3].set_xlabel("M")
    axis[3].set_ylabel("A")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.stimutation.png".format(prefix)), bbox_inches="tight", dpi=300)

    # Scatter of joint posteriors
    post = pd.read_csv("/home/arendeiro/projects/crop-seq/results/digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.ctrl_vs_ctrl.joint_posteriors.csv")
    stats = pd.read_csv(
        os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.ctrl_vs_ctrl.diff_expr.csv")
    )
    stats.index.name = "gene_name"
    genes = stats[abs(stats["cZ"]) > 1.5].index

    fig, axis = plt.subplots(1, figsize=(6, 6))
    axis.hexbin(
        np.log2(1 + post['CROP-seq_Jurkat_TCR_unstimulated']),
        np.log2(1 + post['CROP-seq_Jurkat_TCR_stimulated']),
        alpha=1, bins='log', gridsize=200)
    axis.scatter(
        np.log2(1 + post['CROP-seq_Jurkat_TCR_unstimulated'].ix[genes]),
        np.log2(1 + post['CROP-seq_Jurkat_TCR_stimulated'].ix[genes]),
        alpha=1, color="red")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.joint_posterior.scatter.png".format(prefix)), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.joint_posterior.scatter.svg".format(prefix)), bbox_inches="tight")

    # MA of joint posteriors
    fig, axis = plt.subplots(1, figsize=(6, 6))
    axis.hexbin(
        np.log2(post['CROP-seq_Jurkat_TCR_unstimulated'] * post['CROP-seq_Jurkat_TCR_stimulated']) / 2.,
        stats["Z"].ix[post.index],
        alpha=1, bins='log', gridsize=200)
    axis.scatter(
        np.log2(post['CROP-seq_Jurkat_TCR_unstimulated'].ix[genes] * post['CROP-seq_Jurkat_TCR_stimulated'].ix[genes]) / 2.,
        stats["Z"].ix[genes], color="red", alpha=1)
    axis.set_title("MA plot for {}".format(experiment))
    axis.set_xlabel("M")
    axis.set_ylabel("A")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.joint_posterior.maplot.png".format(prefix)), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.joint_posterior.maplot.svg".format(prefix)), bbox_inches="tight")


def plot_deg_heatmaps(df, assignment, stats, prefix=""):
    """
    """
    from sklearn.manifold import MDS, LocallyLinearEmbedding, Isomap, SpectralEmbedding, TSNE
    from sklearn.decomposition import PCA

    import sys
    sys.setrecursionlimit(5000)

    # remove opposite gRNA library
    if "TCR" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Wnt")]
    elif "WNT" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Tcr")]

    # Add gRNA gene assignment to matrix
    df = df.T
    cell_names = df.index.str.lstrip("un|st")
    cell_ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df["cell_ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in cell_ass], index=df.index).astype("category")
    df["cell_sti"] = pd.Series([x[:2] for x in df.index], index=df.index).astype("category")

    stats = pd.read_csv(
        os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.ctrl_vs_ctrl.diff_expr.csv")
    )
    # remove artitificial chromosomes
    stats = stats[
        (~stats.index.str.startswith("CTRL")) &
        (~stats.index.str.contains("library_"))]

    # fold_pivot = pd.read_csv(os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.diff_expr.csv".format(condition, gene)), index_col=0)

    # get top 500 differential expressed (if significant) in top 2000 cells (by coverage)
    if "WNT" in prefix:
        df2 = df[stats[(stats["q_value"] < 0.05)].sort_values("q_value").head(500).index.tolist() + ["cell_ass", "cell_sti"]]  # filter by p-value
    else:
        df2 = df[stats[(abs(stats["cZ"]) > 1.5)].index.tolist() + ["cell_ass", "cell_sti"]]  # filter by corrected posterior
    # save stats of diff genes
    stats.ix[df2.index.tolist()].to_csv(os.path.join(results_dir, "differential_expression.{}.differential_genes.csv".format(prefix)), index=True)

    # get top 2000 cells from each condition
    df2['sum'] = df2.sum(1)
    big_cells = df2.groupby(['cell_sti'])['sum'].nlargest(1500).index.get_level_values(1)
    df3 = df2.ix[big_cells].drop("sum", axis=1)

    g = sns.clustermap(
        df3._get_numeric_data().T,
        z_score=0,
        # cmap="GnBu",
        robust=True,
        col_colors=get_grna_colors_new(df3, assignment),
        row_colors=get_foldchange_colors_new(df3, stats),
        metric='correlation',
        xticklabels=False, yticklabels=False,
        figsize=(15, 7))
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.clustermap.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.clustermap.svg".format(prefix)), bbox_inches="tight")

    # only assigned cells
    g = sns.clustermap(
        df3[df3["cell_ass"] != "Unassigned"]._get_numeric_data().T,
        z_score=0,
        # cmap="GnBu",
        robust=True,
        col_colors=get_grna_colors_new(df3[df3["cell_ass"] != "Unassigned"], assignment),
        row_colors=get_foldchange_colors_new(df3[df3["cell_ass"] != "Unassigned"], stats),
        metric='correlation',
        xticklabels=False, yticklabels=False,
        figsize=(15, 7))
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.clustermap.only_assigned.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.clustermap.only_assigned.svg".format(prefix)), bbox_inches="tight")

    # sort by stimulation order and gRNA order
    df4 = df3.sort_values(["cell_sti", "cell_ass"])

    g = sns.clustermap(
        df4._get_numeric_data().T,
        z_score=0,
        # cmap="GnBu",
        robust=True,
        col_colors=get_grna_colors_new(df4, assignment),
        row_colors=get_foldchange_colors_new(df4, stats),
        metric='correlation',
        row_cluster=True, col_cluster=False,
        xticklabels=False, yticklabels=True,
        figsize=(15, 7))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sortedcondition_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sortedcondition_heatmap.svg".format(prefix)), bbox_inches="tight")

    # sorted the same but only assigned cells
    z_score = lambda x: (x - x.mean()) / x.std()

    g = sns.clustermap(
        df4[df4["cell_ass"] != "Unassigned"]._get_numeric_data().T.apply(z_score, axis=1).apply(z_score, axis=0),
        # z_score=0,
        robust=True,
        col_colors=get_grna_colors_new(df4[df4["cell_ass"] != "Unassigned"], assignment),
        row_colors=get_foldchange_colors_new(df4[df4["cell_ass"] != "Unassigned"], stats),
        metric='correlation',
        row_cluster=True, col_cluster=False,
        xticklabels=False, yticklabels=True,
        figsize=(15, 7))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sortedconditiongRNA_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sortedconditiongRNA_heatmap.svg".format(prefix)), bbox_inches="tight")

    #

    # Group by stimulation / gRNA, get mean expression

    # use all cells for this (go back to df2)
    df5 = df2.drop("sum", axis=1).groupby(['cell_sti', 'cell_ass']).mean().T

    # report number of cells per condition & gRNA
    fig, axis = plt.subplots(1, figsize=(6, 12))
    c = df2.drop("sum", axis=1).groupby(['cell_sti', 'cell_ass']).apply(len).sort_values(ascending=False)
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
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.clustered.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.clustered.svg".format(prefix)), bbox_inches="tight")

    # sort by stimulation order and gRNA order
    g = sns.clustermap(
        df5.T.sort_index().T, z_score=0,
        col_colors=get_group_colors(df5.T.sort_index().T, assignment),
        row_colors=get_foldchange_colors(df5, stats),
        row_cluster=True, col_cluster=False,
        xticklabels=True, yticklabels=True,
        metric='correlation',
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.sortedconditiongRNA_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.sortedconditiongRNA_heatmap.svg".format(prefix)), bbox_inches="tight")

    #

    # Dimentionality reduction methods
    # try several
    methods = [PCA, LocallyLinearEmbedding, Isomap, SpectralEmbedding, TSNE]

    for name, matrix in [("groups", df5), ("cells", df3.T)]:
        for method in methods:
            print(name, method.__name__)
            model = method()

            if name == "cells":
                color = get_grna_colors(matrix.T, assignment)
                m = matrix.T[matrix.dtypes == np.float64].T
            else:
                color = get_group_colors(matrix.T, assignment)
                m = matrix.T[matrix.dtypes == np.float64]

            fit = model.fit_transform(m)

            # plot
            if method.__name__ == "PCA":
                pcs = 3
            else:
                pcs = 1

            fig, axis = plt.subplots(2, pcs, sharex=False, sharey=False, figsize=(8, 10))
            if method.__name__ != "PCA":
                axis = [[x] for x in axis]

            for i, variable in enumerate(["condition", "gene"]):
                for pc in range(pcs):
                    axis[i][pc].scatter(fit[:, pc], fit[:, pc + 1], color=color[i], alpha=0.75 if name == "groups" else 0.1)
                    axis[i][pc].set_xticklabels([])
                    axis[i][pc].set_yticklabels([])
                    if method.__name__ == "PCA":
                        axis[i][pc].set_xlabel("PC %i" % pc)
                        axis[i][pc].set_ylabel("PC %i" % (pc + 1))
                # axis[i][pc].legend(
                #     handles=[mpatches.Patch(color=v, label=k) for k, v in color_mapping.items()],
                #     ncol=2 if feature == "organ" else 1,
                #     loc='center left',
                #     bbox_to_anchor=(1, 0.5))
            fig.savefig(os.path.join(results_dir, "differential_expression.{}.{}.{}.png".format(prefix, name, method.__name__)), bbox_inches="tight", dpi=300)


def enrich_signature(stats, prefix=""):
    """
    """
    # load stats of diff genes
    stats = pd.read_csv(
        os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.ctrl_vs_ctrl.diff_expr.csv")
    )
    stats.index.name = "gene_name"
    # remove artitificial chromosomes
    stats = stats[
        (~stats.index.str.startswith("CTRL")) &
        (~stats.index.str.contains("library_"))]

    for d, name in [(degs, "_all_genes"), (degs[degs["Z"] > 0], "_up_genes"), (degs[degs["Z"] < 0], "_down_genes")]:

        enr = enrichr(d.reset_index())
        enr.to_csv(os.path.join(results_dir, "differential_expression.{}.enrichr.csv".format(prefix + name)), index=False, encoding="utf8")

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

    # Fig 1c
    fig, axis = plt.subplots()
    sns.distplot(abs(stats["Z"]), kde=False, bins=100, ax=axis)
    # axis.set_yscale("log")
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.posterior_distribution.svg".format(prefix + "_all_genes")), bbox_inches="tight")


def threshold_enrich_signature(stats, prefix=""):
    """
    """
    # load stats of diff genes
    stats = pd.read_csv(
        os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.ctrl_vs_ctrl.diff_expr.csv")
    )
    stats.index.name = "gene_name"
    # remove artitificial chromosomes
    stats = stats[
        (~stats.index.str.startswith("CTRL")) &
        (~stats.index.str.contains("library_"))]

    deg_counts = list()

    for t in np.arange(0.5, 5.0, 0.1):
        sig = stats[abs(stats["cZ"]) > t]
        for d, name in [
                #  (sig, "_all_genes"),
                (sig[sig["Z"] > 0], "_up_genes"), (sig[sig["Z"] < 0], "_down_genes")]:
            deg_counts.append([t, name, d.shape[0]])
            if d.shape[0] > 0:
                enr = enrichr(d.reset_index())
                enr.to_csv(os.path.join(results_dir, "enrichr", "differential_expression.{}.{}.enrichr.csv".format(prefix + name, t)), index=False, encoding="utf8")

    deg_counts = pd.DataFrame(deg_counts)
    deg_counts.columns = ["threshold", "direction", "gene_count"]

    # plot counts of degs
    fig, axis = plt.subplots(len(deg_counts["direction"].unique()), figsize=(6, 8), sharex=True)
    for i, d in enumerate(deg_counts["direction"].unique()):
        p = deg_counts[deg_counts["direction"] == d]
        axis[i].plot(p['threshold'], p['gene_count'], label=d)
        axis[i].set_title(d)
    axis[1].set_ylabel("Number of DEGs")
    axis[2].set_xlabel("threshold")
    fig.savefig(os.path.join(results_dir, "enrichr", "gene_number.svg"), bbox_inches="tight")

    # collect
    enrichments = pd.DataFrame()
    for t in np.arange(0.5, 5.0, 0.1):
        sig = stats[abs(stats["cZ"]) > t]
        for d, name in [(sig[sig["Z"] > 0], "_up_genes"), (sig[sig["Z"] < 0], "_down_genes")]:
            try:
                enr = pd.read_csv(os.path.join(results_dir, "enrichr", "differential_expression.{}.{}.enrichr.csv".format(prefix + name, t)))
            except:
                continue
            enr["threshold"] = t
            enr["name"] = name
            enr["id"] = str(t) + "_" + name
            enrichments = enrichments.append(enr)

    # plot
    for gene_set_library in enrichments["gene_set_library"].unique():
        p = enrichments[enrichments["gene_set_library"] == gene_set_library]

        for metric in ["p_value", "z_score", "combined_score"]:
            enrichments_pivot = pd.pivot_table(p, index="description", columns="id", values=metric).fillna(0)

            # plot
            g = sns.clustermap(
                enrichments_pivot.ix[enrichments_pivot.sum(1).sort_values().tail(50).index],
                col_cluster=False,
                figsize=(15, 15))
            for item in g.ax_heatmap.get_yticklabels():
                item.set_rotation(0)
            for item in g.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
            g.fig.savefig(os.path.join(results_dir, "enrichr", "enrichments.{}.{}.svg".format(gene_set_library, metric)), bbox_inches="tight")


def assign_cells_to_signature(stats, df, assignment, prefix=""):
    """
    """
    # get top 500 differential expressed (if significant) in top 2000 cells (by coverage)
    if "WNT" in prefix:
        stats = stats[(stats["q_value"] < 0.05)]  # filter by p-value
    else:
        stats = stats[(abs(stats["cZ"]) > 1.5)]  # filter by corrected posterior

    # remove opposite gRNA library
    if "TCR" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Wnt")]
    elif "WNT" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Tcr")]

    group_means = pd.read_csv(
        os.path.join(results_dir, "differential_expression.{}.group_means.csv".format(prefix)),
        index_col=0,
        header=[0, 1], skipinitialspace=True, tupleize_cols=True)
    group_means.columns = pd.MultiIndex.from_tuples(group_means.columns)

    # Get stimulation signature
    # 1. get mean expression of each group in signature gens
    x1 = group_means["st", "CTRL"].ix[stats.index]
    x2 = group_means["un", "CTRL"].ix[stats.index]

    # 2. get signature matrix
    # here, bounds are set to (-20, 20) so that the extremes of the signature represent -20% or 120% of the signature
    # this is done because the extreme values (0 and 100%) values correspond to the mean value within each group,
    # meaning that some cells are expected to surpass those values.
    sign = generate_signature_matrix(np.vstack([x1, x2]).T, n=101, bounds=(-20, 20))

    # 3. get signature value of each cell
    sigs = list()
    for i, cell in enumerate(df.columns):
        if i % 100 == 0:
            print(i)
        sigs.append(best_signature_matrix(array=df.ix[stats.index][cell], matrix=sign))

    # Aggregate by condition/gene
    df2 = df.copy().T
    # use all cells for this
    cell_names = df2.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df2["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df2.index).astype("category")
    df2["sti"] = pd.Series([x[:2] for x in df2.index], index=df2.index).astype("category")

    # add signatures
    df2['signature'] = sigs

    df2 = df2[~df2["ass"].isin(["Unassigned", "Essential"])]

    # Save clinical trait signatures for all samples
    df2[["ass", "sti", "signature"]].to_csv(os.path.join(results_dir, "signatures.all_cells.{}scde.csv".format(prefix)))
    df2 = pd.read_csv(os.path.join(results_dir, "signatures.all_cells.{}scde.csv".format(prefix)), index_col=0)

    # plot distibution
    sigs_mean = df2.groupby(['sti', 'ass'])['signature'].median().sort_values().reset_index().set_index(["sti", "ass"])
    sigs_mean["n_cells"] = df2.groupby(['sti', 'ass']).apply(len)

    # Filter out genes with less than the 5th percentile of cells
    sigs_mean = sigs_mean[sigs_mean["n_cells"] >= 10]

    fig, axis = plt.subplots(1, figsize=(10, 8))
    sns.stripplot(x=sigs_mean['signature'], y=sigs_mean.index, orient="horiz", ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.strip.svg".format(prefix)), bbox_inches="tight")

    # Given the known condition of each group, what is the deviation from that?
    # plot as rank of mean
    fig, axis = plt.subplots(1, figsize=(10, 8))
    axis.scatter(sigs_mean['signature'].rank(ascending=False), sigs_mean['signature'])
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.rank.svg".format(prefix)), bbox_inches="tight")

    # plot as violinplots (values per cell)
    p = df2.reset_index().sort_values(['signature'])
    p = p[(p["ass"] != "Essential") & (p["ass"] != "Unassigned")]

    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(
        x="sti", y="signature", hue="ass",
        data=p,
        cut=0,
        order=["st", "un"],
        hue_order=sigs_mean.reset_index().sort_values(["signature"], ascending=False)['ass'].drop_duplicates(),
        ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.violinplot.sorted_st.svg".format(prefix)), bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(
        x="sti", y="signature", hue="ass",
        data=p,
        cut=0,
        order=["st", "un"],
        hue_order=sigs_mean.reset_index().sort_values(["signature"])['ass'].drop_duplicates(),
        ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.violinplot.sorted_un.svg".format(prefix)), bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(x="ass", y="signature", hue="sti", cut=0, data=df2, ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.mean_group_signature.violinplot2.svg".format(prefix)), bbox_inches="tight")

    # Make heatmap sorted by median signature position per knockout

    # Group by stimulation / gRNA, get mean expression

    # use all cells for this
    # diffs = pd.read_csv(os.path.join(results_dir, "differential_expression.{}.differential_genes.csv".format(prefix)), index_col=0)
    df4 = df.copy().T[stats.index].dropna()

    cell_names = df4.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df4["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df4.index).astype("category")
    df4["sti"] = pd.Series([x[:2] for x in df4.index], index=df4.index).astype("category")

    df4 = df4[~df4["ass"].isin(["Unassigned", "Essential"])]

    df5 = df4.groupby(['sti', 'ass']).mean()
    df5 = df5.drop(["Unassigned", "Essential"], level=1)

    # Annotate with median signature and number of cells and save
    df5["signature"] = sigs_mean['signature']
    df5["n_cells"] = df4.groupby(['sti', 'ass']).apply(len).sort_values(ascending=False)
    df5 = df5[df5["n_cells"] >= 10]
    df5.to_csv(os.path.join(results_dir, "differential_expression.{}.group_means.signature_cells_annotated.csv".format(prefix)), index=True)
    df5 = pd.read_csv(os.path.join(results_dir, "differential_expression.{}.group_means.signature_cells_annotated.csv".format(prefix)), index_col=[0, 1], skipinitialspace=True)

    # cluster
    g = sns.clustermap(
        df5.drop(['signature', 'n_cells'], axis=1).T,
        z_score=0,
        col_colors=get_group_colors(df5.drop(['signature', 'n_cells'], axis=1).T, assignment),
        row_colors=get_foldchange_colors_new(df5.drop(['signature', 'n_cells'], axis=1), stats),
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.clustered.png".format(prefix)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.clustered.svg".format(prefix)), bbox_inches="tight")

    # cluster groups, sort genes
    gene_order = (df5.drop(['signature', 'n_cells'], axis=1).ix["st"].median() / df5.drop(['signature', 'n_cells'], axis=1).ix["un"].median()).sort_values(ascending=False).index

    g = sns.clustermap(
        df5[gene_order].T.dropna(), z_score=0,
        col_colors=get_group_colors(df5[gene_order].T, assignment),
        row_colors=get_foldchange_colors_new(df5[gene_order], stats),
        metric='correlation',
        row_cluster=False, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.sorted_genes_clustered_groups.png".format(prefix)), bbox_inches="tight", dpi=300)

    # sort by signature
    p = df5.sort_values("signature").T.drop(['signature', 'n_cells']).ix[gene_order]
    p = p[[i for i, x in enumerate(p.columns.get_level_values(1)) if x not in ["Essential", "Unassigned"]]]

    clust = sns.clustermap(
        p, z_score=0,
        col_colors=get_group_colors(p, assignment),
        row_colors=get_foldchange_colors_new(p.T, stats),
        metric='correlation',
        row_cluster=False, col_cluster=False,
        xticklabels=True, yticklabels=True,
        figsize=(8.62948158106742, 6))
    for item in clust.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in clust.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    clust.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.sorted_signature.png".format(prefix)), bbox_inches="tight", dpi=300)
    clust.fig.savefig(os.path.join(results_dir, "differential_expression.{}.group_means.sorted_signature.svg".format(prefix)), bbox_inches="tight")

    # Strip plot of number of cells
    p = df5.sort_values("signature").T.drop(['signature'])
    p = p[[i for i, x in enumerate(p.columns.get_level_values(1)) if x not in ["Essential", "Unassigned"]]]

    fig, axis = plt.subplots(1, figsize=(8, 8))
    axis.scatter([1] * p.shape[1], range(p.shape[1]), s=p.ix['n_cells'])
    [axis.text(1.0005, i, s=str(int(x))) for i, x in enumerate(p.ix['n_cells'].values)]
    [axis.text(1.002, i, s=str(x)) for i, x in enumerate(p.columns.get_level_values(1))]
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "signatures.all_cells.{}.cells_per_group.bubbles.svg".format(prefix)), bbox_inches="tight")

    # Single-cell matrix sorted in same way as above
    p_all = df4.reset_index().set_index(["sti", "ass"])
    # order gRNAs
    p_all = p_all.ix[p.columns]
    # order genes
    p_all = p_all[[t.get_text() for t in clust.ax_heatmap.get_yticklabels()] + ["index"]]

    # down
    gs = stats[stats["Z"] < 0].sort_values("Z").index.tolist()
    # up
    gs += stats[stats["Z"] > 0].sort_values("Z", ascending=False).index.tolist()
    p_all = p_all[gs + ["index"]]
    # p_all = p_all[diffs.sort_values("fold_change").index.tolist() + ["index"]]

    # for cmap in ["YlGn"]:  # , "YlOrBr", "GnBu", "Greys_r", "Oranges", "ocean_r"]:
    g = sns.clustermap(
        p_all.drop("index", axis=1),
        # standard_scale=0,
        z_score=1,
        # cmap=cmap,
        # vmin=0.25,
        # col_colors=get_foldchange_colors(p_all.set_index("index").T, stats),
        row_colors=get_grna_colors(p_all.set_index("index").T, assignment),
        metric='correlation',
        row_cluster=False, col_cluster=False,
        xticklabels=True, yticklabels=False,
        figsize=(6, 8.62948158106742))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    #g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.assigned_cells.sorted_signature.{}.png".format(prefix, cmap)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.assigned_cells.sorted_signature.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.assigned_cells.sorted_signature.svg".format(prefix)), bbox_inches="tight")

    #

    # Calculate deviation from CTRL for each stimulation
    s = (sigs_mean.ix['st'] - sigs_mean.ix['st', "CTRL"])['signature'].sort_values(ascending=False)
    u = (sigs_mean.ix['un'] - sigs_mean.ix['un', "CTRL"])['signature'].sort_values(ascending=False)
    sl = np.log2(sigs_mean.ix['st'] / sigs_mean.ix['st', "CTRL"])['signature'].sort_values(ascending=False)
    ul = np.log2(sigs_mean.ix['un'] / sigs_mean.ix['un', "CTRL"])['signature'].sort_values(ascending=False)

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


def explore_knockouts(df, assignment, prefix=""):
    """
    """
    df = df.T

    # Annotate cells with stimulus and gRNA assignment
    cell_names = df.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df.index).astype("category")
    df["sti"] = pd.Series([x[:2] for x in df.index], index=df.index).astype("category")

    df2 = df[~df["ass"].isin(["Unassigned", "Essential"])]

    # get condition/gene mean expression for every gene
    group_means = df2.groupby(["sti", "ass"]).mean().dropna().T
    group_means.to_csv(os.path.join(results_dir, "knockout_combination.mean_expression.{}.csv".format(prefix)))

    # filter for groups with more than n cells
    c = df2.groupby(["sti", "ass"]).apply(len)
    group_means = group_means[c[c >= np.percentile(c, 5)].index]

    p = group_means.ix[[x for x in group_means.columns.levels[1].tolist() if x in group_means.index]]
    g = sns.clustermap(
        p,
        robust=True,
        col_colors=get_group_colors(p, assignment),
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.svg".format(prefix)), bbox_inches="tight")

    g = sns.clustermap(
        p,
        z_score=0,
        robust=True,
        col_colors=get_group_colors(p, assignment),
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.z_score.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.z_score.svg".format(prefix)), bbox_inches="tight")

    # Compare difference between conditions
    group_differences = (group_means["st"] - group_means["un"])
    group_differences = group_differences[group_differences.isnull().all()[~group_differences.isnull().all()].index]

    p = group_differences.ix[[x for x in group_differences.columns.tolist() if x in group_differences.index]]
    g = sns.clustermap(
        p,
        robust=True,
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.stimulation_difference.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.stimulation_difference.svg".format(prefix)), bbox_inches="tight")

    g = sns.clustermap(
        p,
        z_score=0,
        robust=True,
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.stimulation_difference.z_score.png".format(prefix)), bbox_inches="tight", dpi=300)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.stimulation_difference.z_score.svg".format(prefix)), bbox_inches="tight")

    # Swarmplots of GATA3, RUNX1, ETS1, EGR1 expression in the two conditions
    fig, axis = plt.subplots(2, figsize=(14, 8))
    sns.violinplot(x="sti", y="GATA3", data=df2, hue="ass")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.GATA3_expression.boxplot.svg".format(prefix)), bbox_inches="tight")
    fig, axis = plt.subplots(2, figsize=(14, 8))
    sns.violinplot(x="sti", y="ETS1", data=df2, hue="ass")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.ETS1_expression.boxplot.svg".format(prefix)), bbox_inches="tight")
    fig, axis = plt.subplots(2, figsize=(14, 8))
    sns.violinplot(x="sti", y="RUNX1", data=df2, hue="ass")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.RUNX1_expression.boxplot.svg".format(prefix)), bbox_inches="tight")
    # Also IL2
    fig, axis = plt.subplots(2, figsize=(14, 8))
    sns.violinplot(x="sti", y="IL2", data=df2, hue="ass")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.own_knockouts.IL2_expression.boxplot.svg".format(prefix)), bbox_inches="tight")

    #

    # All differential genes
    # load all files, select top N significant, append to list
    n = 50
    genes = list()
    for gene in assignment["group"].unique():
        try:
            t = pd.read_csv(os.path.join(results_dir, "differential_expression.{}_onevsctrl.{}.stimutation.csv".format(experiment, gene)))
        except:
            print("Gene {}".format(gene))
            continue
        genes += t.ix[abs(t[t["q_value"] < 0.05]["fold_change"]).sort_values().tail(n).index]["gene"].tolist()
    genes = pd.Series(genes).drop_duplicates()

    p = group_means.ix[genes]
    g = sns.clustermap(
        p,
        robust=True,
        col_colors=get_group_colors(p, assignment),
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.degs.png".format(prefix)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.degs.svg".format(prefix)), bbox_inches="tight")

    g = sns.clustermap(
        p,
        z_score=0,
        robust=True,
        col_colors=get_group_colors(p, assignment),
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.degs.z_score.png".format(prefix)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.degs.z_score.svg".format(prefix)), bbox_inches="tight")

    g = sns.clustermap(
        p[p.sum(0).sort_values().index],
        z_score=0,
        robust=True,
        col_colors=get_group_colors(p, assignment),
        metric='correlation',
        row_cluster=True, col_cluster=False,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.degs.sorted.png".format(prefix)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.mean_expression.degs.z_score.svg".format(prefix)), bbox_inches="tight")

    # load each files, select top X significant, append to list
    enrichments = pd.DataFrame()
    for gene in assignment["group"].unique():
        try:
            t = pd.read_csv(os.path.join(results_dir, "differential_expression.{}_onevsctrl.{}.stimutation.csv".format(experiment, gene)))
        except:
            print("Gene {}".format(gene))
            continue
        print(gene)
        d = t[t["q_value"] < 0.05].rename(columns={"gene": "gene_name"})
        if d.shape[0] < 4:
            continue
        enr = enrichr(d)
        enr["knockout_gene"] = gene

        enrichments = enrichments.append(enr)
    enrichments.to_csv(os.path.join(results_dir, "knockout_combination.{}.degs.enrichr.csv".format(prefix)), index=False, encoding="utf8")
    enrichments = pd.read_csv(os.path.join(results_dir, "knockout_combination.{}.degs.enrichr.csv".format(prefix)))

    for gene_set_library in enrichments["gene_set_library"].unique():

        d = enrichments[enrichments["gene_set_library"] == gene_set_library]

        d_ = pd.pivot_table(d, index="description", columns="knockout_gene", values="combined_score")

        # select terms enriched in most knockouts
        g = sns.clustermap(
            # d_.ix[d_.mean(1).sort_values().tail(35).index].fillna(0),
            d_.ix[(np.nanstd(d_, 1) / d_.sum(1)).sort_values().tail(35).index].fillna(0),
            # z_score=0,
            robust=True,
            row_cluster=True, col_cluster=True,
            xticklabels=True, yticklabels=True,
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
        g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.degs.enrichr.{}.png".format(prefix, gene_set_library)), bbox_inches="tight", dpi=300)
        # g.fig.savefig(os.path.join(results_dir, "knockout_combination.{}.degs.enrichr.csv".format(prefix)), bbox_inches="tight")


def signature_scde(assignment, df):
    # Cluster groups on genes
    df = df.T
    # Annotate cells with stimulus and gRNA assignment
    cell_names = df.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df.index).astype("category")
    df["sti"] = pd.Series([x[:2] for x in df.index], index=df.index).astype("category")

    # read in diff
    degs = pd.read_csv(
        os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.ctrl_vs_ctrl.diff_expr.csv")
    )
    fold_pivot = pd.read_csv(os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.diff_expr.csv".format(condition, gene)), index_col=0)

    # Filter out gRNAs
    degs = degs[~degs.index.str.contains("library|CTRL")]

    # get top N differential genes
    genes = degs[abs(degs["cZ"]) > 1.5].index.tolist()

    df2 = df[df["ass"] == "CTRL"]
    pos = df2[df2["sti"] == "st"].drop(['ass', 'sti'], axis=1).mean(0)
    neg = df2[df2["sti"] == "un"].drop(['ass', 'sti'], axis=1).mean(0)

    # Plot scatter
    fig, axis = plt.subplots(2)
    axis[0].scatter(  # all genes
        neg.ix[neg.index[~neg.index.isin(genes)]],
        pos.ix[pos.index[~pos.index.isin(genes)]],
        s=5,
        color="grey",
        alpha=0.1
    )
    axis[0].scatter(
        neg.ix[genes],
        pos.ix[genes],
        s=10,
        color="orange",
        alpha=0.5
    )
    # Plot MA
    axis[1].scatter(  # all genes
        (neg.ix[neg.index[~neg.index.isin(genes)]] * pos.ix[pos.index[~pos.index.isin(genes)]]) / 2.,
        degs.loc[degs.index[~degs.index.isin(genes)], "cZ"],
        s=5,
        color="grey",
        alpha=0.1
    )
    axis[1].scatter(
        (neg.ix[genes] * pos.ix[genes]) / 2.,
        degs.loc[genes, "cZ"],
        s=10,
        color="orange",
        alpha=0.5
    )

    # Cluster on fold-changes
    g = sns.clustermap(
        fold_pivot.ix[genes].dropna(),
        # robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.fold_change.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)

    # Cluster on original expression
    g = sns.clustermap(
        df.ix[genes].dropna(),
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.expression.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)
    g = sns.clustermap(
        fold_pivot.ix[genes].corr(),
        metric="euclidean",
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.expression.clustering.correlation.png".format(prefix)), bbox_inches="tight", dpi=300)

    # Cluster groups on genes
    df = df.T
    # Annotate cells with stimulus and gRNA assignment
    cell_names = df.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df.index).astype("category")
    df["sti"] = pd.Series([x[:2] for x in df.index], index=df.index).astype("category")

    df2 = df[~df["ass"].isin(["Unassigned", "Essential"])]

    # get condition/gene mean expression for every gene
    group_means = df2.groupby(["sti", "ass"]).mean().dropna().T

    # cluster mean gene expression
    g = sns.clustermap(
        group_means.ix[genes].dropna(),
        z_score=0,
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)
    # correlation
    g = sns.clustermap(
        group_means.ix[genes].dropna().corr(),
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.clustering.correlation.png".format(prefix)), bbox_inches="tight", dpi=300)

    # cluster
    g = sns.clustermap(
        group_means["un"].ix[genes].dropna(),
        z_score=0,
        metric="euclidean",
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.unstimulated.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)
    g = sns.clustermap(
        group_means["st"].ix[genes].dropna(),
        z_score=0,
        metric="euclidean",
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.stimulated.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)


def gather_scde(assignment, N=250):
    """
    """
    fold = pd.DataFrame()

    # remove opposite gRNA library
    if "TCR" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Wnt")]
    elif "WNT" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Tcr")]

    for condition in ["unstimulated"]:
        for gene in assignment[assignment['experiment_group'] == experiment]["group"].unique()[2:-7]:
            if gene in ["Essential", "CTRL"]:
                continue

            # Read up deg
            print(gene)
            try:
                degs = pd.read_csv(
                    os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.CROP-seq_Jurkat_TCR_{}.{}.diff_expr.csv".format(condition, gene))
                )
            except:
                print("Skipping {} {}".format(condition, gene))
                continue
            degs["condition"] = condition
            degs["gene"] = gene
            fold = fold.append(degs)

            enr = enrichr(degs.reset_index().rename(columns={"index": "gene_name"}).head(N))
            enr.to_csv(os.path.join(
                results_dir,
                "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.{}.{}.enrichr.diff_expr.csv".format(condition, gene)), index=False, encoding="utf8")

    # Filter out gRNAs
    fold = fold[~fold.index.str.contains("library")]

    # get top N differential genes
    fold["absZ"] = abs(fold["Z"])
    g = fold.groupby(['condition', 'gene'])['absZ'].nlargest(150)
    genes = g.index.get_level_values(2).unique()

    fold["id"] = fold["condition"] + fold["gene"]

    # create pivot table of fold-changes
    fold_pivot = pd.pivot_table(fold.drop_duplicates().reset_index(), index="index", columns="id", values="Z").dropna()
    fold_pivot.to_csv(os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.diff_expr.csv".format(condition, gene)), index=True)

    # Cluster on fold-changes
    g = sns.clustermap(
        fold_pivot.ix[genes].dropna(),
        # robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        # metric="correlation",
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.fold_change.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)

    # Cluster on original expression
    g = sns.clustermap(
        df.ix[genes].dropna()[np.random.choice(df.columns, 2000)],
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.expression.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)
    g = sns.clustermap(
        fold_pivot.ix[genes].corr(),
        # metric="euclidean",
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.expression.clustering.correlation.png".format(prefix)), bbox_inches="tight", dpi=300)

    # Cluster groups on genes
    df = df.T
    # Annotate cells with stimulus and gRNA assignment
    cell_names = df.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df.index).astype("category")
    df["sti"] = pd.Series([x[:2] for x in df.index], index=df.index).astype("category")

    df2 = df[~df["ass"].isin(["Unassigned", "Essential"])]

    # get condition/gene mean expression for every gene
    group_means = df2.groupby(["sti", "ass"]).mean().dropna().T

    # cluster mean gene expression
    g = sns.clustermap(
        group_means.ix[genes].dropna(),
        z_score=0,
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)
    # correlation
    g = sns.clustermap(
        group_means.ix[genes].dropna().corr(),
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.clustering.correlation.png".format(prefix)), bbox_inches="tight", dpi=300)

    # cluster
    g = sns.clustermap(
        group_means["un"].ix[genes].dropna(),
        z_score=0,
        metric="euclidean",
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.unstimulated.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)
    g = sns.clustermap(
        group_means["st"].ix[genes].dropna(),
        z_score=0,
        metric="euclidean",
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "knockout_signatures.{}.scde.group_expression.stimulated.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)

    # Cluster enrichments
    enrichments = pd.DataFrame()
    for condition in ["stimulated", "unstimulated"]:
        for gene in assignment[assignment['experiment_group'] == experiment]["group"].unique():
            if gene in ["Essential", "CTRL"]:
                continue
            # Read up enrichments
            try:
                enr = pd.read_csv(os.path.join(
                    results_dir,
                    "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.{}.{}.enrichr.diff_expr.csv".format(condition, gene)), encoding="utf8")
            except:
                print("Skipping {} {}".format(condition, gene))
                continue
            enr["condition"] = condition
            enr["gene"] = gene
            enrichments = enrichments.append(enr)
    enrichments["id"] = enrichments["gene"] + enrichments["condition"]
    enrichments.to_csv(os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.enrichr.csv"), index=False, encoding="utf8")
    enrichments = pd.read_csv(os.path.join(results_dir, "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.enrichr.csv"), encoding="utf8")

    for gene_set_library in enrichments["gene_set_library"].unique():

        d = enrichments[enrichments["gene_set_library"] == gene_set_library]

        d_ = pd.pivot_table(d, index="description", columns="id", values="combined_score")

        # cluster
        g = sns.clustermap(
            # terms enriched in most knockouts
            d_.ix[d_.mean(1).sort_values().tail(100).index].fillna(0),
            # most varible terms
            # d_.ix[(np.nanstd(d_, 1) / d_.sum(1)).sort_values().tail(100).index].fillna(0),
            # z_score=0,
            robust=True,
            row_cluster=True, col_cluster=True,
            xticklabels=True, yticklabels=True,
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
        g.fig.savefig(os.path.join(
            results_dir,
            "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.enrichr.{}.svg".format(gene_set_library)), bbox_inches="tight")
        g.fig.savefig(os.path.join(
            results_dir,
            "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.enrichr.{}.pdf".format(gene_set_library)), bbox_inches="tight")

        d_ = pd.pivot_table(d, index="description", columns="id", values="p_value")

        # cluster
        g = sns.clustermap(
            # terms enriched in most knockouts
            d_.ix[d_.mean(1).sort_values().head(100).index].fillna(1),
            # most varible terms
            # d_.ix[(np.nanstd(d_, 1) / d_.sum(1)).sort_values().head(100).index].fillna(1),
            # z_score=0,
            robust=True,
            row_cluster=True, col_cluster=True,
            xticklabels=True, yticklabels=True,
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
        g.fig.savefig(os.path.join(
            results_dir,
            "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.enrichr.{}.p_value.svg".format(gene_set_library)), bbox_inches="tight")
        g.fig.savefig(os.path.join(
            results_dir,
            "digital_expression.500genes.CROP-seq_Jurkat_TCR.scde.enrichr.{}.p_value.pdf".format(gene_set_library)), bbox_inches="tight")


def plot_genes():
    """
    """
    genes = {
        "Cytokines/effector molecules": ["IL2", "TNF", "LTA", "LTB", "GZMB", "IL3", "IL6"],
        "Surface": ["CXCL8", "IL2RA", "CD40LG", "CXCR4", "CD69", "CD83", "TNFRSF4"],
        "maybe chemokines": ["IL8", "CCL4", "CCL3", "CCL2", "CCL1"],
        "other": ["DUSP2", "DUSP5", "VGF", "CSF2", "BIRC3"],
    }
    # matrices
    # raw counts
    exp = pd.read_pickle(os.path.join(results_dir, "digital_expression.{}genes.{}.pickle".format(n_genes, experiment)))

    # raw group counts
    exp2 = exp.to_dense().ix[[y for x in genes.values() for y in x]].dropna().T

    cell_names = exp2.index.str.lstrip("sti:|unsti:")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    exp2['ass'] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=exp2.index).astype("category")
    exp2["sti"] = pd.Series([x[:2] for x in exp2.index], index=exp2.index).astype("category")
    exp_groups = exp2.groupby(['sti', 'ass']).mean().T
    # Filter out genes with less than 10 cells and save
    c = exp2.groupby(['sti', 'ass']).apply(len)
    exp_groups = exp_groups[c[c >= 10].index]
    exp_groups.to_csv(os.path.join(results_dir, "raw_expression.{}.group_means.csv".format(prefix)), index=True)

    # norm counts
    norm = matrix_norm.to_dense()
    # norm group means
    norm2 = norm.ix[[y for x in genes.values() for y in x]].dropna().T

    cell_names = norm2.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    norm2['ass'] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=norm2.index).astype("category")
    norm2["sti"] = pd.Series([x[:2] for x in norm2.index], index=norm2.index).astype("category")
    norm_groups = norm2.groupby(['sti', 'ass']).mean().T
    # Filter out genes with less than 10 cells and save
    c = norm2.groupby(['sti', 'ass']).apply(len)
    norm_groups = norm_groups[c[c >= 10].index]
    norm_groups.to_csv(os.path.join(results_dir, "norm_expression.{}.group_means.csv".format(prefix)), index=True)

    # filter for groups with more than n cells
    c = df2.groupby(["sti", "ass"]).apply(len)
    group_means = group_means[c[c >= np.percentile(c, 5)].index]

    for m, name in [
            (exp_groups, "exp_groups"), (norm_groups, "norm_groups"),  # groups of cells
            # (exp, "counts"), (df, "norm")  # all cells
            ]:
        # take top N covered cells
        g = sns.clustermap(
            m,
            metric="correlation",
            z_score=0)
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
        g.fig.savefig(os.path.join(
            results_dir,
            "gene_expression.500genes.{}.{}.clustermap.svg".format(prefix, name)), bbox_inches="tight")


def classify_unassigned(df, assignment, prefix=""):
    """
    """
    from sklearn.cross_validation import LeaveOneOut
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_curve

    df = df.T

    # Only on differential genes
    diffs = pd.read_csv(os.path.join(results_dir, "differential_expression.{}.differential_genes.csv".format(prefix)), index_col=0)

    # Annotate cells with stimulus and gRNA assignment
    # remove opposite gRNA library
    if "TCR" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Wnt")]
    elif "WNT" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Tcr")]

    cell_names = df.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    df["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df.index).astype("category")
    df["sti"] = pd.Series([x[:2] for x in df.index], index=df.index).astype("category")

    # Split into assigned and unassigned cells
    train = df[~df["ass"].isin(["Unassigned"])][diffs.index.tolist() + ["ass"]]
    predict = df[df["ass"].isin(["Unassigned"])][diffs.index.tolist() + ["ass"]]

    # Approach 1.
    # multiclass learning

    #

    # Approach 2.
    # one vs all learning

    loo = LeaveOneOut(train.shape[0])
    for i, (train_, test_) in enumerate(loo):
        print(i)

        p = (OneVsRestClassifier(RandomForestClassifier(random_state=0))
            .fit(
                train.ix[train_].drop(["ass"], 1),
                train.ix[train_]["ass"])
            .predict_proba(
                train.ix[test_].drop(["ass"], 1)))
        if i == 0:
            scores = p
        else:
            scores = np.concatenate([scores, p])

    scores = pd.DataFrame(scores, columns=sorted(train['ass'].unique()), index=train.index[range(len(scores))])
    scores["truth"] = train["ass"].ix[range(len(scores))]
    scores.to_csv(os.path.join(results_dir, "classification.all_genes.probabilities.csv"))

    # Plot mean group probability for each label
    fig, axis = plt.subplots(1, figsize=(15, 15))
    sns.heatmap(scores.groupby("truth").mean(), cmap="Greens", square=True, ax=axis)
    fig.savefig(os.path.join(results_dir, "classification.all_genes.group_mean_probabilities.svg"), bbox_inches="tight")

    # Make ROC curve (for some classes only maybe)
    roc_curve(train["ass"].ix[range(scores)], scores)
    roc_curve(train["ass"], scores)


def get_grna_colors(d, assignment):
    from matplotlib.colors import colorConverter
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


def get_grna_colors_new(d, assignment):
    from matplotlib.colors import colorConverter
    # make color dict
    pallete = sns.color_palette("colorblind") * 10
    gRNA_color_dict = dict(zip(assignment["group"].unique(), pallete))
    gRNA_color_dict["CTRL"] = "#228B22"
    gRNA_color_dict["Essential"] = "#B22222"
    gRNA_color_dict["Unassigned"] = "grey"

    # match colors with matrix
    stimulation_colors = ["#B22222" if "un" in x else "#228B22" for x in d["cell_sti"]]
    gRNA_colors = [gRNA_color_dict[x] for x in d["cell_ass"]]

    # convert to rgb
    stimulation_colors = [colorConverter.to_rgb(x) if type(x) is str else x for x in stimulation_colors]
    gRNA_colors = [colorConverter.to_rgb(x) if type(x) is str else x for x in gRNA_colors]

    return [stimulation_colors, gRNA_colors]


def get_foldchange_colors_new(d, s):
    fold_change_colors = ["#B22222" if x > 0 else "#333399" for x in s.ix[d.columns]["Z"]]
    return fold_change_colors


def get_group_colors(d, assignment):
    from matplotlib.colors import colorConverter
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
            "NCI-Nature_2016",
            # "Disease_Perturbations_from_GEO_down",
            # "Disease_Perturbations_from_GEO_up",
            # "TF-LOF_Expression_from_GEO"
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


def Wnt_targets():
    """
    """

    # get Wnt target genes
    wnt_targets = pd.read_csv(os.path.join("metadata", "wnt_qPCR_assay_genes.csv")).squeeze()

    # single-cell clustering
    df = matrix_norm.to_dense()
    df2 = df.ix[wnt_targets].dropna()
    # add noise
    noise = np.random.normal(0, 0.0000001, df2.shape[0] * df2.shape[1]).reshape(df2.shape)
    df3 = df2 + noise
    g = sns.clustermap(
        df3,
        xticklabels=False,
        col_colors=["green" if "st" in x else "red" for x in df.columns],
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.single_cell.clustered.png"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(
        df3[sorted(df3.columns)],
        xticklabels=False,
        col_cluster=False,
        col_colors=["green" if "st" in x else "red" for x in df.columns],
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.single_cell.sorted.png"), dpi=300, bbox_inches="tight")

    # z-scored
    g = sns.clustermap(
        df3,
        xticklabels=False,
        col_colors=["green" if "st" in x else "red" for x in df.columns],
        z_score=0,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.single_cell.clustered.zscore.png"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(
        df3[sorted(df3.columns)],
        xticklabels=False,
        col_cluster=False,
        col_colors=["green" if "st" in x else "red" for x in df.columns],
        z_score=0,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.single_cell.sorted.zscore.png"), dpi=300, bbox_inches="tight")

    # group by target clustering
    if "TCR" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Wnt")]
    elif "WNT" in prefix:
        assignment = assignment[~assignment['assignment'].str.contains("Tcr")]

    cell_names = df.index.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])

    df["ass"] = pd.Series([x if type(x) == str else "Unassigned" for x in ass], index=df.index).astype("category")
    df["sti"] = pd.Series([x[:2] for x in df.index], index=df.index).astype("category")

    df2 = df[~df["ass"].isin(["Unassigned", "Essential"])]

    # get condition/gene mean expression for every gene
    group_means = df2.groupby(["sti", "ass"]).mean().dropna().T

    group_means2 = group_means.ix[wnt_targets].dropna()

    g = sns.clustermap(
        group_means2,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.group_means.clustered.png"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(
        group_means2[sorted(group_means2.columns)],
        col_cluster=False,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.group_means.sorted.png"), dpi=300, bbox_inches="tight")

    # z-scored
    g = sns.clustermap(
        group_means2,
        z_score=0,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.group_means.clustered.zscore.png"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(
        group_means2[sorted(group_means2.columns)],
        col_cluster=False,
        z_score=0,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.group_means.sorted.zscore.png"), dpi=300, bbox_inches="tight")

    # CTRL only
    g = sns.clustermap(
        pd.DataFrame([group_means2['st', 'CTRL'], group_means2['un', 'CTRL']]).T,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.group_means.controls.clustered.png"), dpi=300, bbox_inches="tight")
    g = sns.clustermap(
        pd.DataFrame([group_means2['st', 'CTRL'], group_means2['un', 'CTRL']]).T,
        z_score=0,
        metric="correlation"
    )
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "WNT.target_genes.group_means.controls.clustered.zscore.png"), dpi=300, bbox_inches="tight")


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


prj = Project(os.path.join("metadata", "config.separate.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = results_dir = os.path.join("results")

sample_annotation = prj.sheet.df

# get guide annotation
guide_annotation = os.path.join("metadata", "guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)

n_genes = 500
experiment, rows = prj.sheet.df[prj.sheet.df['sample_name'].str.contains("ND")].groupby(['experiment']).groups.items()[0]


s = os.path.join(results_dir, "differential_expression.{}.stimutation.csv".format(prefix))
stats = pd.read_csv(s, index_col=0)

# get expression
for n_genes in [500]:
    for experiment, rows in prj.sheet.df[prj.sheet.df['sample_name'].str.contains("ND")].groupby(['experiment']):
        print(experiment)

        prefix = experiment + "_stimulation.allcells"

        reads = pd.read_csv(os.path.join(results_dir, "{}.guide_cell_gRNA_assignment.all.csv".format(experiment)))
        scores = pd.read_csv(os.path.join(results_dir, "{}.guide_cell_scores.all.csv".format(experiment)))
        assignment = pd.read_csv(os.path.join(results_dir, "{}.guide_cell_assignment.all.csv".format(experiment)))

        exp = pd.read_hdf(os.path.join(prj.paths.results_dir, "{}.digital_expression.{}genes.hdf5".format(experiment, n_genes)), "exp_matrix")

        # Filter matrix
        cells_per_gene = exp.apply(lambda x: (x >= 3).sum(), axis=1)
        genes_per_cell = exp.apply(lambda x: (x >= 3).sum(), axis=0)
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

        hdf5_file = os.path.join(results_dir, "{}.digital_expression.{}genes.seurat_regressed.h5".format(experiment, n_genes))
        matrix_norm = read_seurat_hdf5(hdf5_file, assignment)

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
            stats = differential_genes(pos, neg, assignment, prefix=prefix)
        else:
            stats = pd.read_csv(s, index_col=0)

        if not os.path.exists(os.path.join(results_dir, "differential_expression.{}.stimutation.png".format(prefix))):
            plot_deg_stats(stats, prefix=prefix)

        if not os.path.exists(os.path.join(results_dir, "differential_expression.{}.group_means.png".format(prefix))):
            plot_deg_heatmaps(matrix_norm.to_dense(), assignment, stats, prefix=prefix)

        if not os.path.join(results_dir, "differential_expression.{}.enrichr.csv".format(prefix + "_all_genes")):
            enrich_signature(stats, prefix=prefix)

        assign_cells_to_signature(stats, matrix_norm.to_dense(), assignment, prefix=prefix)

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

        # Aggregate by condition/gene
        # use all cells for this
        for condition in ["st", "un"]:
            for gene in assignment[assignment['experiment_group'] == experiment]["group"].unique():
                if gene in ["Essential", "CTRL"]:
                    continue

                # select condition
                mc = matrix_norm[matrix_norm.columns[matrix_norm.columns.str.contains(condition)].tolist()].to_dense()

                # select cells from this gene
                cell_names = mc.T.index.str.lstrip("un|st")
                ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
                ass = [x if type(x) == str else "Unassigned" for x in ass]

                indices = [i for i, x in enumerate(ass) if x == gene]
                pos = mc.icol(indices)

                # get all other cells
                indices = [i for i, x in enumerate(ass) if x != gene]
                neg = mc.icol(indices)

                prefix = experiment + "_onevsall.{}".format(gene)

                s = os.path.join(results_dir, "differential_expression.{}.stimutation.csv".format(prefix))
                if not os.path.exists(s):
                    stats = differential_genes(pos, neg, assignment, prefix=prefix)
                # else:
                #     stats = pd.read_csv(s, index_col=0)

        # Approach 2:
        # get MSigDB/GEO signatures on these stimuli
        # use them to quantify the deviation of each group of cell with gRNA targeting the same gene

        #

        # Part 2.
        # Explore the knockouts!

        #

        # A.
        # Groupby condition/gene, get mean expression
        # Plot difference/deviation from CTRL for genes:
        # a) of the relevant pathway or
        # b) any differentially expressed gene from the one vs all or one vs CTRL comparisons.
        # b1) with MannU test
        explore_knockouts(matrix_norm.to_dense(), assignment, prefix=prefix)
        # b2) with scde
        gather_scde()
        explore_knockouts_scde(, assignment

        # Plot the difference between stimulated/unstimulated for same genes

        # B.
        # Get enrichments for all DEGs from the one vs all or one vs CTRL comparisons.
        # Plot knockout vs enrichment
