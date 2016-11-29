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

import sys
sys.setrecursionlimit(10000)


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
            srun -J seurat --mem 256000 -p longq -n 1 -c 24 --pty /cm/shared/apps/R/3.2.3/bin/R

            function(to_norm, output){
                library(Seurat)
                library(dplyr)
                library(Matrix)

                library("rhdf5")
                root_dir = "/home/arendeiro/projects/crop-seq/results"
                hdf_file = "CROP-seq_Jurkat_TCR.digital_expression.500genes.only_assigned.hdf5.gz"
                hdf_block = "exp_matrix"
                output_prefix = "CROP-seq_Jurkat_TCR.digital_expression.500genes.scde"

                root_dir = path.expand(root_dir)  # hdf5 needs absolute paths

                # Load the exp dataset
                raw_counts = h5read(paste(root_dir, hdf_file, sep="/"), hdf_block)

                # VERY IMPORTANT!: remember R is 1-based and the level indexes come from Python!
                condition = as.character(raw_counts$block0_items_level0[raw_counts$block0_items_label0 + 1])
                replicate = as.factor(as.character(raw_counts$block0_items_level1[raw_counts$block0_items_label1 + 1]))
                cell = as.character(raw_counts$block0_items_level2[raw_counts$block0_items_label2 + 1])
                grna = as.character(raw_counts$block0_items_level3[raw_counts$block0_items_label3 + 1])

                # Make dataframe
                M = as.data.frame(t(raw_counts$block0_values))
                M <-apply(M, 2, function(x) {storage.mode(x) <- 'integer'; x})
                names = paste(condition, replicate, cell, grna, sep="|")
                colnames(M) <- names
                rownames(M) <- raw_counts$axis1

                # Initialize the Seurat object with the raw data
                exp <- new("seurat", raw.data = M)

                # Keep all genes expressed in at least 20 cells, keep all cells with >= 500 genes
                exp <- Setup(
                    exp, project="cropseq",
                    min.cells=5, min.genes=500,
                    do.logNormalize=T, total.expr=1e4)

                # Stash cell attributes
                exp <- AddMetaData(exp, condition, "condition")
                exp <- AddMetaData(exp, as.factor(replicate), "replicate")
                exp <- AddMetaData(exp, cell, "cell")
                exp <- AddMetaData(exp, grna, "grna")

                # Calculate the percentage of mitochondrial genes and store it.
                mito.genes <- grep("^MT-", rownames(exp@data), value=T)
                percent.mito <- colSums(expm1(exp@data[mito.genes, ])) / colSums(expm1(exp@data))
                exp <- AddMetaData(exp, percent.mito, "percent.mito")
                # Calculate the percentage of ribosomal genes and store it.
                ribo.genes <- grep("^RP", rownames(exp@data), value=T)
                percent.ribo <- colSums(expm1(exp@data[ribo.genes, ])) / colSums(expm1(exp@data))
                exp <- AddMetaData(exp, percent.ribo, "percent.ribo")

                # # Plot QC stats
                # pdf(paste(root_dir, "seurat.qc_plots.pdf", sep="/"))
                # # VlnPlot(exp, c("replicate", "nGene", "nUMI", "percent.mito"))
                # # GenePlot(exp, "IL2RA", "MALAT1", cex.use = 1)
                # par(mfrow = c(2, 2))
                # GenePlot(exp, "nUMI", "nGene")
                # GenePlot(exp, "nUMI", "percent.mito")
                # GenePlot(exp, "nUMI", "percent.ribo")
                # GenePlot(exp, "percent.ribo", "percent.mito")
                # dev.off()

                # Filter out cells with more than 10% of mitochondrial transcripts
                exp_subset <- SubsetData(
                    exp, subset.name="percent.mito",
                    accept.high=0.10)  # quantile(percent.mito, 0.95) # 95th percentile of mitochondrial transcripts

                # Regress out
                exp_subset_regressed <- RegressOut(exp_subset, latent.vars=c("nUMI", "percent.mito"))

                # Write as hdf5
                output_name = "CROP-seq_Jurkat_TCR.digital_expression.500genes.only_assigned.seurat_regressed.hdf5.gz"
                output = paste(root_dir, output_name, sep="/")
                h5createFile(output)
                h5createGroup(output, "seurat_matrix")
                h5write(as.matrix(exp_subset_regressed@scale.data), file=output, "seurat_matrix/matrix")
                h5write(colnames(exp_subset_regressed@scale.data), file=output, "seurat_matrix/columns")
                h5write(rownames(exp_subset_regressed@scale.data), file=output, "seurat_matrix/rows")
                H5close()

                return(as.matrix(exp_subset_regressed@data))


                # continue to explore
                exp_subset_regressed <- MeanVarPlot(exp_subset_regressed ,fxn.x=expMean, fxn.y=logVarDivMean, x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5, do.contour=F, do.plot=F)

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


def read_seurat_hdf5(hdf5_file):
    import h5py
    with h5py.File(hdf5_file, 'r') as handle:
        cols = handle.get("seurat_matrix/columns").value
        rows = handle.get("seurat_matrix/rows").value
        df = handle.get("seurat_matrix/matrix").value
    seurat_matrix = pd.DataFrame(df, index=cols, columns=rows).T

    # add info as multiindex columns
    condition = map(lambda x: x[0], seurat_matrix.columns.str.split("|"))
    replicate = map(lambda x: x[1], seurat_matrix.columns.str.split("|"))
    cell = map(lambda x: x[2], seurat_matrix.columns.str.split("|"))
    grna = map(lambda x: x[3], seurat_matrix.columns.str.split("|"))
    gene = map(lambda x: x[1] if len(x) > 1 else x[0][:4], pd.Series(grna).str.split("_"))
    seurat_matrix.columns = pd.MultiIndex.from_arrays([condition, replicate, cell, grna, gene], names=['condition', 'replicate', 'cell', 'grna', 'gene'])

    return seurat_matrix


def unsupervised(df, prefix=""):
    # Inspect
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE, MDS, LocallyLinearEmbedding, SpectralEmbedding, Isomap
    from matplotlib import cm

    methods = ["PCA", "TSNE", "LocallyLinearEmbedding", "SpectralEmbedding", "Isomap", "MDS"]

    fig, axis = plt.subplots(len(df.columns.levels), len(methods), figsize=(4 * len(df.columns.levels), 4 * len(methods)))
    for i, method in enumerate(methods):
        fitted = eval(method)().fit_transform(df)

        for j, level in enumerate(df.columns.names):
            print(method, level)

            # color mapping
            integer_map = dict([(val, cm.viridis(i)) for i, val in enumerate(set(df.columns.get_level_values(level)))])
            colors = [integer_map[x] for x in df.columns.get_level_values(level)]

            axis[i, j].scatter(fitted[:, 0], fitted[:, 1], color=colors, alpha=0.1)
    fig.savefig(os.path.join(results_dir, "clustering.{}.png".format(prefix)), bbox_inches="tight", dpi=300)


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


def stimulation_signature(assignment, df, experiment="CROP-seq_Jurkat_TCR", n_genes=500, cond1="stimulated", cond2="unstimulated"):
    def get_level_colors(index):
        from matplotlib.colors import colorConverter
        pallete = sns.color_palette("colorblind") * int(1e6)

        colors = list()

        if hasattr(index, "levels"):
            for level in index.levels:
                color_dict = dict(zip(level, pallete))
                level_colors = [color_dict[x] for x in index.get_level_values(level.name)]
                colors.append(level_colors)
        else:
            color_dict = dict(zip(set(index), pallete))
            index_colors = [color_dict[x] for x in index]
            colors.append(index_colors)

        return colors

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

        res = dict()
        for i in range(matrix.shape[1]):
            res[i] = pearsonr(array, matrix[:, i])

        cors = {i: x[0] for i, x in res.items()}
        p_values = {i: x[1] for i, x in res.items()}
        return cors, p_values
        # return cors.values().index(max(cors.values()))  # index
        # (
        # cors.values().index(max(cors.values())),  # index
        # max(cors.values())  # highest correlation value


    # read in diff
    degs = pd.read_csv(
        os.path.join(results_dir, "{}.digital_expression.{}genes.scde.diff_expr.csv".format(experiment, n_genes)), index_col=0
    )
    degs[cond1] = np.log2(degs[cond1 + "_posterior"] + 1)
    degs[cond2] = np.log2(degs[cond2 + "_posterior"] + 1)

    # Filter out gRNAs
    degs = degs[~degs.index.str.contains("library|CTRL")]

    # get top N differential genes
    de_genes = degs[abs(degs["cZ"]) > 2].index.tolist()
    degs["sig"] = (abs(degs["cZ"]) > 2)
    o_genes = df.index[~df.index.isin(genes)]
    degs["A"] = np.log2(degs[cond1 + "_posterior"] * degs[cond2 + "_posterior"]) / 2.

    # df2 = df[df.columns[df.columns.get_level_values('gene') == "CTRL"]]
    # pos = df2[df.columns.get_level_values('condition') == "stimulated"].mean(0)
    # neg = df2[df.columns.get_level_values('condition') == "unstimulated"].mean(0)

    # Plot scatter of group posteriors
    fig, axis = plt.subplots(1, figsize=(4, 4 * 1))
    axis.scatter(  # all genes
        degs[cond1].ix[o_genes],
        degs[cond2].ix[o_genes],
        s=5,
        color="grey",
        alpha=0.1
    )
    axis.scatter(
        degs[cond1].ix[de_genes],
        degs[cond2].ix[de_genes],
        s=10,
        color="orange",
        alpha=0.5
    )
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.scatter.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.scatter.svg".format(experiment, n_genes)), bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(4, 4 * 1))
    # Plot MA
    axis.scatter(  # all genes
        degs.loc[o_genes, "A"],
        degs.loc[o_genes, "Z"],
        s=5,
        color="grey",
        alpha=0.1
    )
    axis.scatter(
        degs.loc[de_genes, "A"],
        degs.loc[de_genes, "Z"],
        s=10,
        color="orange",
        alpha=0.5
    )
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.maplot.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.maplot.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Cluster all cells on DE genes
    g = sns.clustermap(
        df.ix[de_genes],
        metric="correlation",
        # robust=True
        vmin=-4, vmax=4,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.all_cells.clustering.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # Correlate all cells on DE genes
    g = sns.clustermap(
        df.ix[de_genes].corr(),
        metric="correlation",
        # robust=True
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df.columns),
        row_colors=get_level_colors(df.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.all_cells.clustering.correlation.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # Cluster CTRL cells on DE genes
    g = sns.clustermap(
        df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].ix[de_genes],
        # robust=True
        metric="correlation",
        vmin=-4, vmax=4,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.control_cells.clustering.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # Correlate CTRL cells on DE genes
    g = sns.clustermap(
        df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].ix[de_genes].corr(),
        metric="correlation",
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].columns),
        row_colors=get_level_colors(df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.control_cells.clustering.correlation.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # Cluster groups of targeted gRNAs/genes
    # get condition/gene mean expression for every gene
    df_grna_means = df.T.groupby(level=["condition", "grna"]).mean().T
    df_gene_means = df.T.groupby(level=["condition", "gene"]).mean().T

    # grna level
    # cluster mean gene expression
    g = sns.clustermap(
        df_grna_means.ix[de_genes],
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        col_colors=get_level_colors(df_grna_means.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.grna_means.clustering.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # order
    g = sns.clustermap(
        df_grna_means.ix[de_genes],
        robust=True,
        row_cluster=True, col_cluster=False,
        yticklabels=False, xticklabels=True,
        col_colors=get_level_colors(df_grna_means.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.grna_means.ordered_heatmap.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # correlation
    g = sns.clustermap(
        df_grna_means.ix[de_genes].corr(),
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        col_colors=get_level_colors(df_grna_means.columns),
        row_colors=get_level_colors(df_grna_means.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.grna_means.clustering.correlation.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # gene level
    # cluster mean gene expression
    g = sns.clustermap(
        df_gene_means.ix[de_genes],
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        col_colors=get_level_colors(df_grna_means.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.gene_means.clustering.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # order
    g = sns.clustermap(
        df_gene_means.ix[de_genes],
        robust=True,
        row_cluster=True, col_cluster=False,
        yticklabels=False, xticklabels=True,
        col_colors=get_level_colors(df_gene_means.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.gene_means.ordered_heatmap.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # correlation
    g = sns.clustermap(
        df_gene_means.ix[de_genes].corr(),
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        col_colors=get_level_colors(df_grna_means.columns),
        row_colors=get_level_colors(df_grna_means.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.gene_means.clustering.correlation.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    #

    #

    # Signature-based cell assignemnt

    # Get stimulation signature
    # 1. get mean expression of each group in signature gens
    df_grna_means = df.T.groupby(level=["condition", "grna"]).mean().T
    df_gene_means = df.T.groupby(level=["condition", "gene"]).mean().T
    x1 = df_gene_means[df_gene_means.columns[df_gene_means.columns.get_level_values("condition") == cond1]].mean(axis=1).ix[de_genes]
    x2 = df_gene_means[df_gene_means.columns[df_gene_means.columns.get_level_values("condition") == cond2]].mean(axis=1).ix[de_genes]

    # 2. get signature matrix
    # here, bounds are set to (-20, 20) so that the extremes of the signature represent -20% or 120% of the signature
    # this is done because the extreme values (0 and 100%) values correspond to the mean value within each group,
    # meaning that some cells are expected to surpass those values.
    sign_mat = generate_signature_matrix(np.vstack([x1, x2]).T, n=101, bounds=(-20, 20))

    # 3. get signature value of each cell
    cors = list()
    p_values = list()
    for i, cell in enumerate(df.columns):
        if i % 100 == 0:
            print(i)
        cor, p_value = best_signature_matrix(array=df.ix[de_genes][cell], matrix=sign_mat)
        cors.append(cor)
        p_values.append(p_value)

    cors = pd.DataFrame(cors, index=df.columns)
    p_values = pd.DataFrame(p_values, index=df.columns)

    # 4. get background of signature correlations/positions
    n_permutations = 10
    for i in range(n_permutations):
        background_matrix = df.ix[de_genes].copy().values
        np.random.shuffle(background_matrix)
        c = list()
        p = list()
        for i in range(background_matrix.shape[1]):
            if i % 100 == 0:
                print(i)
            cor, p_value = best_signature_matrix(array=background_matrix[:, i], matrix=sign_mat)
            c.append(cor)
            p.append(p_value)

        if i == 0:
            random_cors = pd.DataFrame(c)
            random_p_values = pd.DataFrame(p)
        else:
            random_cors = random_cors.append(pd.DataFrame(c))
            random_p_values = random_p_values.append(pd.DataFrame(p))

    # 5. investigate signatures
    sigs = cors.apply(lambda x: np.argmax(x), axis=1)
    sigs.name = "signature"
    sigs_p = p_values.apply(lambda x: np.argmax(x), axis=1)
    sigs_r = random_cors.apply(lambda x: np.argmax(x), axis=1)
    sigs_r.name = "signature"
    sigs_rp = random_p_values.apply(lambda x: np.argmax(x), axis=1)

    # get "uncorrelated" fraction
    res = 1 - cors.max(axis=1)
    res.name = "residual"
    sigs = pd.DataFrame(sigs).join(res)
    res_r = 1 - sigs_r
    res_r.name = "residual"
    sigs_r = pd.DataFrame(sigs_r).join(res_r)

    # save
    sigs.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation.csv".format(experiment, n_genes)))

    # all cells together vs random
    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4))
    sns.distplot(sigs, ax=axis[0])
    sns.distplot(p_values.apply(lambda x: np.argmin(x), axis=1), ax=axis[0])
    sns.distplot(random_cors.apply(lambda x: np.argmax(x), axis=1), ax=axis[1])
    sns.distplot(random_p_values.apply(lambda x: np.argmin(x), axis=1), ax=axis[1])
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_pvalue.distributions.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # per gene vs random
    g = sns.FacetGrid(
        data=cors.apply(lambda x: np.argmax(x), axis=1).reset_index(),
        row='condition', col='gene')
    g.map(sns.violinplot, "signature")
    g.map(sns.stripplot, "signature")
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.per_gene.correlation.distributions.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # signature assignemnt vs residual
    g = sns.FacetGrid(
        data=sigs.reset_index(),
        col='condition')
    g.map(plt.scatter, "signature", "residual", alpha=0.2)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.correlation.signature_vs_residual.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # residual vs read coverage per cell
    g = sns.FacetGrid(
        data=sigs.reset_index(),
        col='condition')
    g.map(plt.scatter, "signature", "residual", alpha=0.2)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.correlation.signature_vs_residual.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    g = sns.FacetGrid(
        data=sigs.reset_index(),
        row='condition', col='gene')
    g.map(plt.scatter, "signature", "residual", alpha=0.2)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.per_gene.correlation.signature_vs_residual.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)


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


prj = Project(os.path.join("metadata", "config.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = results_dir = os.path.join("results")

sample_annotation = prj.sheet.df

# get guide annotation
guide_annotation = os.path.join("metadata", "guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)

# n_genes = 500
# experiment, rows = prj.sheet.df.groupby(['experiment']).groups.items()[1]


# get expression
for n_genes in [500]:
    for experiment, rows in prj.sheet.df.groupby(['experiment']):
        print(experiment)

        prefix = experiment + "_stimulation.allcells"

        assignment = pd.read_csv(os.path.join(results_dir, "{}.guide_cell_assignment.all.csv".format(experiment)))

        # exp = pd.read_hdf(os.path.join(results_dir, "{}.digital_expression.{}genes.hdf5.gz".format(experiment, n_genes)), "exp_matrix", compression="gzip")
        exp_assigned = pd.read_hdf(os.path.join(results_dir, "{}.digital_expression.{}genes.only_assigned.hdf5.gz".format(experiment, n_genes)), "exp_matrix", compression="gzip")

        # Normalize
        # Approach 1:
        # normalize by total
        # matrix_norm = normalize(matrix, experiment=experiment, kind="total")

        # Approach 2:
        # regress out based on total number and MT-genes using Seurat
        matrix_norm = normalize(exp_assigned, kind="seurat")

        hdf5_file = os.path.join(results_dir, "{}.digital_expression.{}genes.only_assigned.seurat_regressed.hdf5.gz".format(experiment, n_genes))
        matrix_norm = read_seurat_hdf5(hdf5_file)

        #

        # Unsupervised

        # Approach 1:
        # apply dimentionality reduction methods/clustering
        # and discover biological axis related with stimulation
        unsupervised(matrix_norm)

        # Approach 2 (just a bit supervised though):
        # get differential genes between conditions from CTRL cells only
        # use this signature to position each cell
        # observe deviation of groups of cells with gRNA targeting the same
        # gene coming from either condition to the mean

        # variant A:

        # get differential genes with scde
        # stats = differential_genes(pos, neg, assignment, prefix=prefix) <- todo: add scde R code wrapped here

        # visualize signature and make signature position assignments per cell/gRNA/gene
        stimulation_signature(assignment, matrix_norm)

        # Part 2.
        # Explore the knockouts!

        #

        # A.
        # Groupby condition/gene, get mean expression
        # Plot difference/deviation from CTRL for genes:
        # a) of the relevant pathway or
        # b) any differentially expressed gene from the one vs all or one vs CTRL comparisons.
        # b1) with MannU test
        # explore_knockouts(matrix_norm.to_dense(), assignment, prefix=prefix)
        # b2) with scde
        # gather_scde()
        # explore_knockouts_scde()

        # Plot the difference between stimulated/unstimulated for same genes

        # B.
        # Get enrichments for all DEGs from the one vs all or one vs CTRL comparisons.
        # Plot knockout vs enrichment
