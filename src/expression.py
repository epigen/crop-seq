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


def differential_genes(pos, neg, df, assignment, prefix=""):
    """
    """
    from scipy.stats import mannwhitneyu
    import multiprocessing
    import parmap
    from statsmodels.sandbox.stats.multicomp import multipletests

    def test(g):
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
                    test,
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

    # get top 500 differential expressed (if significant) in top 2000 cells (by coverage)
    df2 = df.ix[abs(stats["fold_change"]).sort_values().head(500).index.tolist()][df.sum().sort_values().tail(2000).index.tolist()].to_dense()

    # color gRNAs
    # remove opposite library
    if "TCR" in experiment:
        assignment = assignment[~assignment['assignment'].str.contains("Wnt")]
    elif "WNT" in experiment:
        assignment = assignment[~assignment['assignment'].str.contains("Tcr")]
    # make color dict
    pallete = sns.color_palette("colorblind") * 10
    gRNA_color_dict = dict(zip(assignment["group"].unique(), pallete))
    gRNA_color_dict[pd.np.nan] = "grey"

    # match colors with current sort order
    stimulation_colors = ["red" if "un" in x else "green" for x in df2.columns]
    cell_names = df2.columns.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    ass = [x if type(x) == str else pd.np.nan for x in ass]
    gRNA_colors = [gRNA_color_dict[x] for x in ass]

    g = sns.clustermap(
        df2, z_score=0,
        col_colors=[stimulation_colors, gRNA_colors],
        cmap="GnBu",
        xticklabels=False, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.clustermap.png".format(prefix)), bbox_inches="tight", dpi=300)

    # sort order
    e = df2.ix[stats.fold_change.sort_values().index].dropna()
    e = e[e.columns[e.sum() != 0]]
    ee = e.T.sort_values(e.index.tolist(), ascending=False).T

    # match colors with current sort order
    stimulation_colors = ["red" if "un" in x else "green" for x in ee.columns]
    cell_names = ee.columns.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    ass = [x if type(x) == str else pd.np.nan for x in ass]
    gRNA_colors = [gRNA_color_dict[x] for x in ass]

    g = sns.clustermap(
        ee, z_score=0,
        col_colors=[stimulation_colors, gRNA_colors],
        cmap="GnBu",
        row_cluster=False, col_cluster=False,
        xticklabels=False, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sorted_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)

    # sort by stimulation order
    df3 = df2.sort_index(axis=1)
    # match colors with current sort order
    stimulation_colors = ["red" if "un" in x else "green" for x in df3.columns]
    cell_names = df3.columns.str.lstrip("un|st")
    ass = pd.Series([assignment[assignment["cell"] == y]["group"].squeeze() for y in cell_names])
    ass = [x if type(x) == str else pd.np.nan for x in ass]
    gRNA_colors = [gRNA_color_dict[x] for x in ass]

    g = sns.clustermap(
        ee, z_score=0,
        col_colors=[stimulation_colors, gRNA_colors],
        cmap="GnBu",
        row_cluster=True, col_cluster=False,
        xticklabels=False, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.sortedcondition_heatmap.png".format(prefix)), bbox_inches="tight", dpi=300)

    #

    #

    #


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
experiment = "CROP-seq_HEK293T_4_WNT"

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

        differential_genes(pos, neg, matrix_norm, assignment, prefix=experiment + "_stimulation.allcells")

        # variant B:
        # only "CTRL" cells
        experiment_assignment = assignment[assignment["experiment_group"] == experiment]
        ctrl_cells = experiment_assignment[experiment_assignment['group'] == "CTRL"]['cell']
        cell_names = matrix_norm.columns.str.lstrip("un|st")

        selected_cells = pd.Series(matrix_norm.columns).where(cell_names.isin(ctrl_cells)).dropna().tolist()

        ctrl_matrix = matrix_norm[selected_cells]

        pos = ctrl_matrix[ctrl_matrix.columns[ctrl_matrix.columns.str.contains("st")].tolist()].to_dense()
        neg = ctrl_matrix[ctrl_matrix.columns[ctrl_matrix.columns.str.contains("un")].tolist()].to_dense()

        differential_genes(pos, neg, ctrl_matrix, assignment, prefix=experiment + "_stimulation.ctrlcells")

        # Supervised
        # Approach 1:
        # Compare groups of cells with gRNA targeting the same gene with either:
        # a) control cells; b) all cells; - intersect (or max) the differential genes from either,
        # explore signature enrichment

        # Approach 2:
        # get MSigDB/GEO signatures on these stimuli
        # use them to quantify the deviation of each group of cell with gRNA targeting the same gene
