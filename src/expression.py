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
        return np.log2(1 + df.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)

    def seurat_regress(df):
        """
        """
        import rpy2.robjects as robj
        import pandas.rpy.common as com

        run = robj.r("""
            # module unload R
            # module load gcc/6.0.0
            # module load texlive
            # module load R/3.2.3
            # srun -J seurat --mem 256000 -p longq -n 1 -c 24 --pty --x11 /cm/shared/apps/R/3.2.3/bin/R

            function(hdf_file, output_file){
                library(Seurat)
                library(dplyr)
                library(Matrix)

                library("rhdf5")

                # root_dir = "/home/arendeiro/projects/crop-seq/results"
                # hdf_file = "CROP-seq_Jurkat_TCR.digital_expression.500genes.only_assigned.hdf5.gz"
                # hdf_file = paste(root_dir, hdf_file, sep="/")
                # hdf_block = "exp_matrix"
                # root_dir = path.expand(root_dir)  # hdf5 needs absolute paths
                # output_name = "CROP-seq_Jurkat_TCR.digital_expression.500genes.only_assigned.seurat_regressed.hdf5.gz"
                # output_file = paste(root_dir, output_name, sep="/")

                # Load the exp dataset
                raw_counts = h5read(hdf_file, hdf_block)

                # VERY IMPORTANT!: remember R is 1-based and the level indexes come from Python (0-based)!
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
                h5createFile(output_file)
                h5createGroup(output_file, "seurat_matrix")
                h5write(as.matrix(exp_subset_regressed@scale.data), file=output_file, "seurat_matrix/matrix")
                h5write(colnames(exp_subset_regressed@scale.data), file=output_file, "seurat_matrix/columns")
                h5write(rownames(exp_subset_regressed@scale.data), file=output_file, "seurat_matrix/rows")
                H5close()
            }
        """)
        run(df,
            os.path.join(results_dir, "{}.digital_expression.500genes.only_assigned.seurat_regressed.csv".format(experiment)),
            os.path.join(results_dir, "{}.digital_expression.500genes.only_assigned.seurat_regressed.hdf5.gz".format(experiment)))

        return read_seurat_hdf5(os.path.join(results_dir, "{}.digital_expression.500genes.only_assigned.seurat_regressed.hdf5.gz".format(experiment)))

    if kind == "total":
        norm = normalize_by_total(df)
        norm.to_hdf(os.path.join(results_dir, "digital_expression.500genes.{}.log2_tpm.hdf5.gz".format(experiment)), "log2_tpm", compression="gzip")
        return norm
    else:
        regressed = seurat_regress(df).to_sparse(fill_value=0)
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


def unsupervised(df, experiment="", filter_low=True):
    # Inspect
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE, MDS, LocallyLinearEmbedding, SpectralEmbedding, Isomap

    methods = ["PCA", "TSNE", "LocallyLinearEmbedding", "SpectralEmbedding", "Isomap", "MDS"]

    if filter_low:
        df = df[df.sum().sort_values().tail(int(df.shape[1] * 0.25)).index]  # get top 25% covered cells

    level_mapping = dict(zip(df.columns.names, range(len(df.columns.names))))

    # Cells grouped by gene or grna
    for group in ["CTRL_cells", "all_cells"]:
        df2 = df.copy()
        if group == 'CTRL_cells':
            df2 = df2[df2.columns[df2.columns.get_level_values("gene") == "CTRL"]]

        for level in ["gene", "grna"]:
            df_group = df2.T.groupby(level=[level_mapping['condition'], level_mapping[level]]).median().T
            df_group.index = df2.index

            fig, axis = plt.subplots(len(df_group.columns.levels), len(methods), figsize=(3 * len(methods), 3 * len(df_group.columns.levels)))
            for i, method in enumerate(methods):
                model = eval(method)()
                try:
                    fitted = model.fit_transform(df_group.T)
                except:
                    continue

                for j, level in enumerate(df_group.columns.names):
                    if level == "cell":
                        continue
                    print(group, method, level)

                    # color mapping
                    integer_map = dict(zip(df_group.columns.levels[j], sns.color_palette("colorblind") * int(1e5)))
                    colors = [integer_map[x] for x in df_group.columns.get_level_values(level)]

                    axis[j, i].scatter(fitted[:, 0], fitted[:, 1], color=colors, alpha=0.4)

                    if j == len(df_group.columns.levels) - 1:
                        for p in range(fitted.shape[0]):
                            axis[j, i].text(fitted[p, 0], fitted[p, 1], " ".join(df_group.columns[p]), color=colors[p], alpha=0.8)

                    if method == "PCA":
                        fig2, axis2 = plt.subplots(1)
                        axis2.plot(range(1, fitted.shape[0] + 1), (model.explained_variance_ / model.explained_variance_.sum()) * 100, "-o")
                        axis2.set_xlabel("PC")
                        axis2.set_ylabel("% variance explained")
                        fig2.savefig(os.path.join(results_dir, "{}.clustering.{}.grouped_{}.PCA_variance.svg".format(experiment, group, level)), bbox_inches="tight")
            # fig.savefig(os.path.join(results_dir, "{}.clustering.{}.grouped_{}.png".format(experiment, group, level)), bbox_inches="tight", dpi=300)
            fig.savefig(os.path.join(results_dir, "{}.clustering.{}.grouped_{}.svg".format(experiment, group, level)), bbox_inches="tight")

    methods = ["PCA", "TSNE", "LocallyLinearEmbedding", "SpectralEmbedding", "Isomap"]

    for group in ["CTRL_cells", "all_cells"]:
        if group == 'CTRL_cells':
            df2 = df2[df.columns[df.columns.get_level_values("gene") == "CTRL"]]
        else:
            df2 = df

        fig, axis = plt.subplots(len(df2.columns.levels), len(methods) + 1, figsize=(3 * len(df2.columns.levels), 3 * len(methods)))
        for i, method in enumerate(methods):
            model = eval(method)()
            fitted = model.fit_transform(df2.T)

            for j, level in enumerate(df2.columns.names):
                if level == "cell":
                    continue
                print(method, level)

                # color mapping
                integer_map = dict(zip(df2.columns.levels[level_mapping[level]], sns.color_palette("colorblind") * int(1e5)))
                colors = [integer_map[x] for x in df2.columns.get_level_values(level)]

                axis[j, i].scatter(fitted[:, 0], fitted[:, 1], color=colors, alpha=0.4)

                if method == "PCA":
                    fig2, axis2 = plt.subplots(1)
                    axis2.plot(range(1, fitted.shape[0] + 1), (model.explained_variance_ / model.explained_variance_.sum()) * 100, "-o")
                    axis2.set_xlabel("PC")
                    axis2.set_ylabel("% variance explained")
                    fig2.savefig(os.path.join(results_dir, "{}.clustering.{}.PCA_variance.svg".format(experiment, group)), bbox_inches="tight")

                    if level == "condition":
                        fig3, axis3 = plt.subplots(4, 4, figsize=(4 * 3, 4 * 3))
                        for pc in range(4 * 4):
                            axis3.flat[pc].scatter(fitted[:, pc], fitted[:, pc + 1], color=colors, alpha=0.4)
                            axis3.flat[pc].set_xlabel(pc + 1)
                            axis3.flat[pc].set_ylabel(pc + 2)
                        fig3.savefig(os.path.join(results_dir, "{}.clustering.{}.PCA.PCs.svg".format(experiment, group)), bbox_inches="tight")

        fig.savefig(os.path.join(results_dir, "{}.clustering.{}.png".format(experiment, group)), bbox_inches="tight", dpi=300)


def get_level_colors(index):
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


def differential_genes(df, experiment, method="pca"):
    """
    """
    def scde():
        return

    def deseq2(df):
        df.to_csv("count_matrix.csv")
        pd.DataFrame([df.columns, df.columns.str.get(0)], index=["cell", "condition"]).T.to_csv("col_matrix.csv", index=False)
        return

    def pca(df, level="gene", filter_low=True):
        from sklearn.decomposition import PCA
        level_mapping = dict(zip(df.columns.names, range(len(df.columns.names))))

        if filter_low:
            df2 = df[df.sum().sort_values().tail(int(df.shape[1] * 0.25)).index]  # get top 25% covered cells
        else:
            df2 = df

        # Cells grouped by gene
        df_group = df2.T.groupby(level=[level_mapping["condition"], level_mapping[level]]).median().T
        df_group.index = df2.index

        fitted = PCA().fit_transform(df_group)  # for genes

        r = pd.Series(fitted[:, 1], index=df2.index).sort_values()
        g = sns.clustermap(
            df2.ix[r[abs(r) > np.percentile(abs(r), 99)].index],
            metric="correlation",
            z_score=0,
            vmin=-3, vmax=3,
            row_cluster=True, col_cluster=True,
            yticklabels=True, xticklabels=False,
            col_colors=get_level_colors(df2.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.PCA.clustering.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

        r = pd.Series(fitted[:, 1], index=df2.index).sort_values()
        g = sns.clustermap(
            df_group.ix[r[abs(r) > np.percentile(abs(r), 95)].index],
            metric="correlation",
            z_score=0,
            vmin=-3, vmax=3,
            row_cluster=True, col_cluster=True,
            yticklabels=True, xticklabels=True,
            col_colors=get_level_colors(df_group.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.PCA.clustering.{}_level.png".format(experiment, n_genes, level)), bbox_inches="tight", dpi=300)

        return r

    print("Getting differential gene expression for experiment: '{}'".format(experiment))

    if method == "scde":
        scde()
    elif method == "deseq":
        deseq2(df)
    elif method == "pca":
        diff = pca(df)
        diff.to_csv(os.path.join(results_dir, "{}.differential_expression.{}.stimutation.csv".format(experiment, method)), index=True)


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


def stimulation_signature(df, experiment="CROP-seq_Jurkat_TCR", n_genes=500, method="pca", cond1="stimulated", cond2="unstimulated"):

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

    if method == "scde":
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
        pd.Series(de_genes).to_csv("~/degenes.csv", index=False)
        o_genes = df.index[~df.index.isin(de_genes)]
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
        axis.set_xlabel(cond1)
        axis.set_ylabel(cond2)
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
        axis.set_xlabel("A")
        axis.set_ylabel("M {} vs {}".format(cond1, cond2))
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.maplot.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
        fig.savefig(os.path.join(results_dir, "{}.{}genes.scde_dge.maplot.svg".format(experiment, n_genes)), bbox_inches="tight")

    # if using pca genes
    if method == "pca":
        diff = pd.read_csv(os.path.join(results_dir, "{}.differential_expression.{}.stimutation.csv".format(experiment, method)), squeeze=True, index_col=0, header=None, names=["gene_name"])
        de_genes = diff[abs(diff) > np.percentile(abs(diff), 99)].index.tolist()

    # Cluster all cells on DE genes
    g = sns.clustermap(
        df.ix[de_genes],
        metric="correlation",
        z_score=0,
        # robust=True
        vmin=-3, vmax=3,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.{}.all_cells.clustering.png".format(experiment, n_genes, method)), bbox_inches="tight", dpi=300)

    # Correlate all cells on DE genes
    g = sns.clustermap(
        df.ix[de_genes].corr(),
        vmin=0, vmax=1, cmap="BrBG",
        metric="correlation",
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df.columns),
        row_colors=get_level_colors(df.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.{}.all_cells.clustering.correlation.png".format(experiment, n_genes, method)), bbox_inches="tight", dpi=300)

    # Cluster CTRL cells on DE genes
    g = sns.clustermap(
        df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].ix[de_genes],
        metric="correlation",
        z_score=0,
        # robust=True
        vmin=-3, vmax=3,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df.columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.{}.control_cells.clustering.png".format(experiment, n_genes, method)), bbox_inches="tight", dpi=300)

    # Correlate CTRL cells on DE genes
    g = sns.clustermap(
        df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].ix[de_genes].corr(),
        vmin=0, vmax=1, cmap="BrBG",
        metric="correlation",
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        col_colors=get_level_colors(df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].columns),
        row_colors=get_level_colors(df[df.columns[df.columns.get_level_values('gene') == "CTRL"]].columns),
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
        item.set_fontsize(8)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.{}.control_cells.clustering.correlation.png".format(experiment, n_genes, method)), bbox_inches="tight", dpi=300)

    # Cluster groups of targeted gRNAs/genes
    # get condition/gene mean expression for every gene
    for level in ["gene", "grna"]:
        df_grna_means = df.T.groupby(level=["condition", level]).mean().T
        # grna level
        # cluster mean gene expression
        g = sns.clustermap(
            df_grna_means.ix[de_genes],
            z_score=0,
            vmin=-3, vmax=3,
            # robust=True,
            row_cluster=True, col_cluster=True,
            yticklabels=False, xticklabels=True,
            col_colors=get_level_colors(df_grna_means.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.{}.{}_level.clustering.png".format(experiment, n_genes, method, level)), bbox_inches="tight", dpi=300)

        # order
        g = sns.clustermap(
            df_grna_means.ix[de_genes],
            z_score=0,
            vmin=-3, vmax=3,
            # robust=True,
            row_cluster=True, col_cluster=False,
            yticklabels=False, xticklabels=True,
            col_colors=get_level_colors(df_grna_means.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.{}.{}_level.ordered_heatmap.png".format(experiment, n_genes, method, level)), bbox_inches="tight", dpi=300)

        # correlation
        g = sns.clustermap(
            df_grna_means.ix[de_genes].corr(),
            vmin=0, vmax=1, cmap="BrBG",
            row_cluster=True, col_cluster=True,
            yticklabels=True, xticklabels=False,
            col_colors=get_level_colors(df_grna_means.columns),
            row_colors=get_level_colors(df_grna_means.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.{}.{}_level.clustering.correlation.png".format(experiment, n_genes, method, level)), bbox_inches="tight", dpi=300)

    #

    # Signature-based cell assignemnt

    # Get stimulation signature
    # 1. get mean expression of each group in signature gens
    df_grna_means = df.T.groupby(level=["condition", "grna"]).median().T
    df_gene_means = df.T.groupby(level=["condition", "gene"]).median().T
    x1 = df_gene_means[df_gene_means.columns[df_gene_means.columns.get_level_values("condition") == cond1]].median(axis=1).ix[de_genes]
    x1.name = cond1
    x2 = df_gene_means[df_gene_means.columns[df_gene_means.columns.get_level_values("condition") == cond2]].median(axis=1).ix[de_genes]
    x2.name = cond2

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

    cors.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_correlation.csv".format(experiment, n_genes)))
    p_values.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_p_values.csv".format(experiment, n_genes)))

    # visualize
    # clustered
    g = sns.clustermap(
        cors,
        z_score=0,
        row_cluster=True, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(cors.index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_matrix.clustered.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        -np.log10(p_values),
        z_score=0,
        row_cluster=True, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(p_values.index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.pvalue_matrix.clustered.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    # sorted by max
    sigs = cors.apply(lambda x: np.argmax(x), axis=1).sort_values()

    g = sns.clustermap(
        cors.ix[sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(cors.ix[sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_matrix.sorted.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        -np.log10(p_values).ix[sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(p_values.ix[sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.pvalue_matrix.sorted.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    for level in ['grna', 'gene']:
        c = cors.groupby(level=['condition', level]).mean()
        cs = c.apply(lambda x: np.argmax(x), axis=1).sort_values()

        p = p_values.groupby(level=['condition', level])

        g = sns.clustermap(
            c,
            z_score=0,
            row_cluster=True, col_cluster=False,
            yticklabels=True, xticklabels=False,
            row_colors=get_level_colors(c.index))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_matrix.clustered.{}.png".format(experiment, n_genes, level)), dpi=300, bbox_inches="tight")

        g = sns.clustermap(
            c.ix[cs.index],
            z_score=0,
            row_cluster=False, col_cluster=False,
            yticklabels=True, xticklabels=False,
            row_colors=get_level_colors(c.ix[cs.index].index))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_pvalue_matrix.sorted.{}.png".format(experiment, n_genes, level)), dpi=300, bbox_inches="tight")

    # 4. get background of signature correlations/positions
    n_permutations = 100
    for i in range(n_permutations):
        background_matrix = df.ix[de_genes].copy().values.flatten()
        np.random.shuffle(background_matrix)  # shuffle only shuffles in one dimention!
        background_matrix = background_matrix.reshape(df.ix[de_genes].shape)
        c = list()
        p = list()
        for j in range(background_matrix.shape[1]):
            if j % 100 == 0:
                print(i, j)
            cor, p_value = best_signature_matrix(array=background_matrix[:, j], matrix=sign_mat)
            c.append(cor)
            p.append(p_value)

        if i == 0:
            random_cors = pd.DataFrame(c, index=df.columns)
            random_p_values = pd.DataFrame(p, index=df.columns)
        else:
            random_cors = random_cors.append(pd.DataFrame(c, index=df.columns))
            random_p_values = random_p_values.append(pd.DataFrame(p, index=df.columns))

    random_cors.to_csv(os.path.join(results_dir, "{}.{}genes.signature.random_cells.matrix_correlation.csv".format(experiment, n_genes)))
    random_p_values.to_csv(os.path.join(results_dir, "{}.{}genes.signature.random_cells.matrix_p_values.csv".format(experiment, n_genes)))

    # visualize
    # clustered
    g = sns.clustermap(
        random_cors,
        z_score=0,
        row_cluster=True, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(random_cors.index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.correlation_matrix.clustered.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        -np.log10(random_p_values),
        z_score=0,
        row_cluster=True, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(random_p_values.index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.pvalue_matrix.clustered.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    # sorted by max
    random_sigs = random_cors.apply(lambda x: np.argmax(x), axis=1).sort_values()

    g = sns.clustermap(
        random_cors.ix[random_sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(random_cors.ix[random_sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.correlation_matrix.sorted.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        -np.log10(random_p_values).ix[random_sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(random_p_values.ix[random_sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.pvalue_matrix.sorted.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    for level in ['grna', 'gene']:
        c = random_cors.groupby(level=['condition', level]).mean()
        cs = c.apply(lambda x: np.argmax(x), axis=1).sort_values()

        p = random_p_values.groupby(level=['condition', level])

        g = sns.clustermap(
            c,
            z_score=0,
            row_cluster=True, col_cluster=False,
            yticklabels=True, xticklabels=False,
            row_colors=get_level_colors(c.index))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.correlation_matrix.clustered.{}.png".format(experiment, n_genes, level)), dpi=300, bbox_inches="tight")

        g = sns.clustermap(
            c.ix[cs.index],
            z_score=0,
            row_cluster=False, col_cluster=False,
            yticklabels=True, xticklabels=False,
            row_colors=get_level_colors(c.ix[cs.index].index))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.matrix_pvalue_matrix.sorted.{}.png".format(experiment, n_genes, level)), dpi=300, bbox_inches="tight")

    # # Try quadratic programming
    # fits_raw = list()
    # for i, cell in enumerate(df.columns):
    #     if i % 100 == 0:
    #         print(i)
    #     fit = quadratic_signature(array=df.ix[de_genes].dropna()[cell], matrix=pd.DataFrame(x1).join(x2))
    #     fits_raw.append(fit)

    # fits = pd.DataFrame(fits_raw, index=df.columns)
    # fits = fits.join(fits['x'].apply(pd.Series)).drop('x', axis=1).rename(columns={0: cond1, 1: cond2})

    # 5. investigate signatures
    sigs = cors.apply(lambda x: np.argmax(x), axis=1)
    sigs.name = "signature"
    sigs_r = random_cors.apply(lambda x: np.argmax(x), axis=1)
    sigs_r.name = "signature"

    # get "uncorrelated" fraction (noise)
    res = 1 - cors.max(axis=1)
    res.name = "residual"
    res_r = 1 - sigs_r
    res_r.name = "residual"

    # reads per cell
    sigs = pd.merge(sigs.reset_index(), res.reset_index()).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])

    # reads per cell
    s = exp_assigned.sum(axis=0)
    s.name = "reads_per_cell"
    sigs = pd.merge(sigs.reset_index(), s.reset_index()).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])

    # get minimum correlation
    m = cors.min(axis=1)
    m.name = "min_corr"
    sigs = pd.merge(sigs.reset_index(), m.reset_index()).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])

    # get mean correlation
    m = cors.mean(axis=1)
    m.name = "mean_corr"
    sigs = pd.merge(sigs.reset_index(), m.reset_index()).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])

    # # get "smoothness" of ranked cumulative sum
    # def ranked_cumsum(x, rank_threshold=0.25):
    #     # scale from 0 to 1
    #     xs = (x - x.min()) / (x.max() - x.min())
    #     rs = xs.sort_values(ascending=False).cumsum() / xs.sum()
    #     return (rs < rank_threshold).sum()

    # r = cors.apply(ranked_cumsum, axis=1)

    # save
    sigs.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation.csv".format(experiment, n_genes)))
    sigs = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation.csv".format(experiment, n_genes))).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])
    # sigs['replicate'] = sigs['replicate'].astype(str)

    # pairwise variable distribution
    g = sns.pairplot(sigs.reset_index(), vars=sigs.columns, hue="condition")
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_metrics.pairwise_distributions.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # per gene vs random
    g = sns.FacetGrid(
        data=cors.apply(lambda x: np.argmax(x), axis=1).reset_index(),
        row='condition', col='gene')
    g.map(sns.violinplot, "signature")
    g.map(sns.stripplot, "signature")
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.per_gene.correlation.distributions.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # all cells together vs random
    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4))
    sns.distplot(sigs, ax=axis[0])
    sns.distplot(p_values.apply(lambda x: np.argmin(x), axis=1), ax=axis[0])
    sns.distplot(random_cors.apply(lambda x: np.argmax(x), axis=1), ax=axis[1])
    sns.distplot(random_p_values.apply(lambda x: np.argmin(x), axis=1), ax=axis[1])
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_pvalue.distributions.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # annotate KO genes with signature
    sigs_mean = sigs.groupby(level=['condition', 'gene']).mean()
    sigs_mean["n_cells"] = sigs.groupby(level=['condition', 'gene']).apply(len)

    # add distance from CTRL
    for cond in [cond1, cond2]:
        ko = sigs_mean.loc[sigs_mean.index.get_level_values('condition') == cond, "signature"]
        ctrl = sigs_mean.loc[(sigs_mean.index.get_level_values('condition') == cond) & (sigs_mean.index.get_level_values('gene') == "CTRL"), "signature"]
        sigs_mean.loc[sigs_mean.index.get_level_values('condition') == cond, "abs_change"] = ko - ctrl.squeeze()
        sigs_mean.loc[sigs_mean.index.get_level_values('condition') == cond, "log_fold_change"] = np.log2(ko / ctrl.squeeze())

    # save
    sigs_mean = sigs_mean.sort_index(level='condition')
    sigs_mean = sigs_mean.sort_values(['signature'])
    sigs_mean.to_csv(os.path.join(results_dir, "{}.{}genes.signature.group_means.annotated.csv".format(experiment, n_genes)), index=True)

    #

    # Filter out genes with less than n of cells
    sigs_mean = sigs_mean[sigs_mean["n_cells"] >= 10]

    fig, axis = plt.subplots(1, figsize=(10, 8))
    sns.stripplot(x=sigs_mean['signature'], y=sigs_mean.index, orient="horiz", ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.all_cells.mean_group_signature.strip.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Given the known condition of each group, what is the deviation from that?
    # plot as rank of mean
    fig, axis = plt.subplots(1, figsize=(10, 8))
    axis.scatter(sigs_mean['signature'].rank(ascending=False), sigs_mean['signature'])
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.all_cells.mean_group_signature.rank.svg".format(experiment, n_genes)), bbox_inches="tight")

    # plot as violinplots (values per cell)
    p = sigs.reset_index().sort_values(['signature'])

    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(
        x="condition", y="signature", hue="gene",
        data=p,
        cut=0,
        order=[cond1, cond2],
        hue_order=sigs_mean.reset_index().sort_values(["signature"], ascending=False)['gene'].drop_duplicates(),
        ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.all_cells.mean_group_signature.violinplot.sorted_st.svg".format(experiment, n_genes)), bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(
        x="condition", y="signature", hue="gene",
        data=p,
        cut=0,
        order=[cond1, cond2],
        hue_order=sigs_mean.reset_index().sort_values(["signature"])["gene"].drop_duplicates(),
        ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.all_cells.mean_group_signature.violinplot.sorted_un.svg".format(experiment, n_genes)), bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(16, 8))
    sns.violinplot(x="gene", y="signature", hue="condition", cut=0, data=sigs.reset_index(), ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.all_cells.mean_group_signature.violinplot2.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Make heatmap sorted by median signature position per knockout

    # Group by stimulation / gRNA, get mean expression
    mean_df = df.ix[de_genes].T.groupby(level=['condition', 'gene']).mean()

    # cluster
    g = sns.clustermap(
        mean_df.T, z_score=0,
        col_colors=get_level_colors(mean_df.index),
        metric='correlation',
        row_cluster=True, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.gene_group_means.clustered.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression..group_means.clustered.svg".format(experiment, n_genes)), bbox_inches="tight")

    # cluster groups, sort genes
    g = sns.clustermap(
        mean_df.ix[sigs_mean.sort_values("signature").index].T, z_score=0,
        col_colors=get_level_colors(mean_df.ix[sigs_mean.sort_values("signature").index].index),
        metric='correlation',
        row_cluster=False, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.gene_group_means.sorted_genes_clustered_groups.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # sort by signature
    clust = sns.clustermap(
        mean_df.ix[sigs_mean.sort_values("signature").index], z_score=1,
        row_colors=get_level_colors(sigs_mean.sort_values("signature").index),
        metric='correlation',
        row_cluster=False, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(8.62948158106742, 6))
    for item in clust.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in clust.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    clust.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.gene_group_means.sorted_signature.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
    clust.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.gene_group_means.sorted_signature.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Strip plot of number of cells
    p = sigs_mean.sort_values("signature")

    fig, axis = plt.subplots(1, figsize=(8, 8))
    axis.scatter([1] * p.shape[0], range(p.shape[0]), s=p['n_cells'])
    [axis.text(1.0005, i, s=str(int(x))) for i, x in enumerate(p['n_cells'].values)]
    [axis.text(1.01, i, s=str(x)) for i, x in enumerate(p.index.get_level_values('condition'))]
    [axis.text(1.02, i, s=str(x)) for i, x in enumerate(p.index.get_level_values('gene'))]
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.cells_per_group.bubbles.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Stripe with relative change
    fig, axis = plt.subplots(1, figsize=(8, 8))
    sns.heatmap(p[["log_fold_change"]], ax=axis)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.cells_per_group.stripe.svg".format(experiment, n_genes)), bbox_inches="tight")

    # # Single-cell matrix sorted in same way as above
    m = p['signature'].to_dict()

    df2 = df.copy().T
    df2.index = df2.index.droplevel(['cell', 'replicate', 'grna'])
    df2['sortby'] = df2.index.to_series().map(m).tolist()
    df2 = df2.sort_values("sortby").drop('sortby', axis=1).T
    df2.columns = df.columns

    # sort by signature
    g = sns.clustermap(
        df2.ix[de_genes].T, z_score=1,
        row_colors=get_level_colors(df2.columns),
        metric='correlation',
        row_cluster=False, col_cluster=True,
        xticklabels=True, yticklabels=True,
        figsize=(8.62948158106742, 6))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.single_cells.sorted_signature.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.single_cells.sorted_signature.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Barplots of difference compared with CTRL
    fig, axis = plt.subplots(1, figsize=(4, 4))
    sns.barplot(sigs_mean['abs_change'], sigs_mean['abs_change'].index, orient="horiz", order=sigs_mean['abs_change'].index, ax=axis[0])
    axis[0].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_ylabel(cond1)
    axis[1].set_ylabel(cond2)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.mean_group_signature_deviation.ranklog2.barplot.svg".format(experiment, n_genes)), bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis[0].scatter(sigs_mean['log_fold_change'].rank(ascending=False), sigs_mean['log_fold_change'])
    axis[0].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_ylabel(cond1)
    axis[1].set_ylabel(cond2)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.mean_group_signature_deviation.ranklog2.scatter.svg".format(experiment, n_genes)), bbox_inches="tight")


def gather_scde(assignment, N=30, experiment="CROP-seq_Jurkat_TCR", n_genes=500, conds=["stimulated", "unstimulated"]):
    """
    """
    def get_level_colors(index):
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

    scde_output = pd.DataFrame()

    for condition in conds:
        for gene in assignment["gene"].drop_duplicates():
            if gene in ["library", "CTRL"]:
                continue

            # Read up deg
            print(condition, gene)
            try:
                degs = pd.read_csv(
                    os.path.join(results_dir, "{}.digital_expression.{}genes.scde.diff_expr.{}.{}.csv".format(experiment, n_genes, condition, gene))
                )
            except IOError:
                print("Skipping {} {}".format(condition, gene))
                continue
            degs["experiment"] = experiment
            degs["n_genes"] = n_genes
            degs["condition"] = condition
            degs["gene"] = gene
            degs = degs.rename(columns={gene + "_posterior": "gene_posterior"})
            scde_output = scde_output.append(degs.reset_index())

            # enr = enrichr(degs.reset_index().rename(columns={"index": "gene_name"}).head(N))
            # enr.to_csv(os.path.join(
            #     results_dir,
            #     "{}.digital_expression.{}genes.scde.diff_expr.{}.{}.enrichr.csv".format(experiment, n_genes, condition, gene)), index=False, encoding="utf8")

    # Calculate A and if sig
    scde_output["A"] = np.log2(scde_output["gene_posterior"] * scde_output["CTRL_posterior"]) / 2.
    scde_output["sig"] = (abs(scde_output["cZ"]) > 2)

    scde_output = scde_output.set_index("index")

    # Filter out gRNAs
    scde_output = scde_output[~scde_output.index.str.contains("library|_gene")]

    scde_output.to_csv(os.path.join(
        results_dir,
        "{}.digital_expression.{}genes.scde.knockouts.all_conditions_genes.csv".format(experiment, n_genes)), index=True)

    scde_output = pd.read_csv(os.path.join(
        results_dir,
        "{}.digital_expression.{}genes.scde.knockouts.all_conditions_genes.csv".format(experiment, n_genes)), index_col=0)

    # get top N differential genes
    scde_output["absZ"] = abs(scde_output["Z"])
    g = scde_output.groupby(['condition', 'gene'])['absZ'].nlargest(N)
    genes = g.index.get_level_values("index").unique()

    scde_output["id"] = scde_output["condition"] + scde_output["gene"]

    # create pivot table of fold-changes
    scde_pivot = pd.pivot_table(scde_output.drop_duplicates().reset_index(), index="index", columns="id", values="Z").dropna()
    scde_pivot.to_csv(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.csv".format(experiment, n_genes, N)), index=True)
    scde_pivot = pd.read_csv(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.csv".format(experiment, n_genes, N)), index_col=0)

    # Cluster on fold-changes
    g = sns.clustermap(
        scde_pivot.ix[genes].dropna(),
        # robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        metric="correlation",
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.fold_change.clustering.png".format(experiment, n_genes, N)), bbox_inches="tight", dpi=300)

    # Cluster on original expression
    g = sns.clustermap(
        matrix_norm.ix[genes],
        robust=True,
        metric="correlation",
        col_colors=get_level_colors(matrix_norm.columns),
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.expression.clustering.png".format(experiment, n_genes, N)), bbox_inches="tight", dpi=300)

    g = sns.clustermap(
        scde_pivot.ix[genes].corr(),
        metric="correlation",
        col_colors=get_level_colors(matrix_norm.columns),
        row_colors=get_level_colors(matrix_norm.columns),
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.expression.clustering.correlation.png".format(experiment, n_genes, N)), bbox_inches="tight", dpi=300)

    # Cluster groups on genes
    mean_matrix_norm = matrix_norm.T.groupby(level=['condition', 'gene']).mean().T

    # gene vs gene matrix
    self_genes = assignment['gene'].drop_duplicates().sort_values()
    self_genes = self_genes[~self_genes.isin(["CTRL"])]

    g = sns.clustermap(
        mean_matrix_norm.ix[self_genes],
        metric="correlation",
        col_colors=get_level_colors(mean_matrix_norm.columns),
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.self_genes.group_expression.clustering.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    for cond in conds:
        p = mean_matrix_norm[mean_matrix_norm.columns[mean_matrix_norm.columns.get_level_values('condition') == cond]]
        g = sns.clustermap(
            p.ix[self_genes],
            z_score=0,
            col_colors=get_level_colors(p.columns),
            row_cluster=False, col_cluster=False,
            yticklabels=True, xticklabels=True,
            figsize=(15, 15))
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.self_genes.group_expression.{}.heatmap.png".format(experiment, n_genes, condition)), bbox_inches="tight", dpi=300)

    # cluster mean gene expression
    g = sns.clustermap(
        mean_matrix_norm.ix[genes],
        metric="correlation",
        col_colors=get_level_colors(matrix_norm.columns),
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.group_expression.clustering.png".format(experiment, n_genes, N)), bbox_inches="tight", dpi=300)

    # correlation
    g = sns.clustermap(
        mean_matrix_norm.ix[genes].corr(),
        metric="correlation",
        col_colors=get_level_colors(matrix_norm.columns),
        row_colors=get_level_colors(matrix_norm.columns),
        row_cluster=True, col_cluster=True,
        yticklabels=True, xticklabels=False,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.group_expression.clustering.correlation.png".format(experiment, n_genes, N)), bbox_inches="tight", dpi=300)

    # cluster
    m1 = mean_matrix_norm.T.ix[mean_matrix_norm.T.index[mean_matrix_norm.T.index.get_level_values('condition') == conds[0]]].T
    g = sns.clustermap(
        m1.ix[genes],
        col_colors=get_level_colors(m1.columns),
        metric="correlation",
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.group_expression.{}.clustering.png".format(experiment, n_genes, N, conds[0])), bbox_inches="tight", dpi=300)

    m2 = mean_matrix_norm.T.ix[mean_matrix_norm.T.index[mean_matrix_norm.T.index.get_level_values('condition') == conds[0]]].T
    g = sns.clustermap(
        m2.ix[genes],
        col_colors=get_level_colors(m2.columns),
        metric="correlation",
        robust=True,
        row_cluster=True, col_cluster=True,
        yticklabels=False, xticklabels=True,
        figsize=(15, 15))
    for item in g.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    g.fig.savefig(os.path.join(results_dir, "{}.digital_expression.{}genes.scde.knockouts.top{}.group_expression.{}.clustering.png".format(experiment, n_genes, N, conds[1])), bbox_inches="tight", dpi=300)

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


def plot_genes(exp):
    """
    """
    genes = {
        "Cytokines/effector molecules": ["IL2", "TNF", "LTA", "LTB", "GZMB", "IL3", "IL6"],
        "Surface": ["CXCL8", "IL2RA", "CD40LG", "CXCR4", "CD69", "CD83", "TNFRSF4"],
        "maybe chemokines": ["IL8", "CCL4", "CCL3", "CCL2", "CCL1"],
        "other": ["DUSP2", "DUSP5", "VGF", "CSF2", "BIRC3"],
    }
    # matrices
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


def inspect_bulk():
    """
    """
    def normalize_quantiles_r(dataframe):
        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()

        robjects.r('require("preprocessCore")')
        normq = robjects.r('normalize.quantiles')

        dataframe_norm = pd.DataFrame(
            np.array(normq(dataframe)),
            index=dataframe.index,
            columns=dataframe.columns
        )

        return dataframe_norm

    # read gene expression matrix
    bitseq = pd.read_csv(os.path.join("results", "{}.count_matrix.gene_level.csv".format(experiment)), index_col=[0], header=range(4))
    esat = pd.read_csv(os.path.join("results", "{}.ESAT_count_matrix.csv".format(experiment)), index_col=[0], header=range(4)).fillna(0)

    bitseq_quant = normalize_quantiles_r(bitseq)
    esat_quant = normalize_quantiles_r(esat)

    bitseq = np.log2(1 + bitseq.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)
    esat = np.log2(1 + esat.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)

    # bitseq_seurat = seurat_regress(bitseq)
    # esat_seurat = seurat_regress(esat)
    bitseq_seurat = pd.read_csv(os.path.join("results", "{}.count_matrix.gene_level.seurat_regressed.csv".format(experiment)), index_col=[0])
    bitseq_seurat.columns = bitseq.columns
    esat_seurat = pd.read_csv(os.path.join("results", "{}.ESAT_count_matrix.seurat_regressed.csv".format(experiment)), index_col=[0])
    esat_seurat.columns = esat.columns

    # read in diff genes
    degs = pd.read_csv(
        os.path.join(results_dir, "{}.digital_expression.{}genes.scde.diff_expr.csv".format(experiment, n_genes)), index_col=0
    )
    degs = degs[~degs.index.str.contains("library|CTRL")]
    de_genes = degs[abs(degs["cZ"]) > 2].index.tolist()

    # read in diff genes from Bulk
    degs = pd.read_csv(
        os.path.join(results_dir, "bulk", "deseq_esat_CTRL_samples_stimulation_signature.diff.csv"), index_col=0
    )
    de_genes_bulk = degs[degs["padj"] < 0.05].index.tolist()

    quant_types = [
        ("bitseq", bitseq), ("esat", esat), ("bitseq_quant", bitseq_quant),
        ("esat_quant", esat_quant), ("bitseq_seurat", bitseq_seurat), ("esat_seurat", esat_seurat)]

    for quant_type, exp_matrix in quant_types:
        print(quant_type)
        # exp_matrix = exp_matrix.ix[matrix_norm.index].dropna()
        # matrix_norm = matrix_norm.ix[exp_matrix.index].dropna()
        # de_genes = [x for x in de_genes if x in exp_matrix.index]

        # Boxplots of expression
        fig, axis = plt.subplots(1)
        sns.boxplot(data=pd.melt(exp_matrix), x="grna", y="value", hue="condition", ax=axis)
        fig.savefig(os.path.join("results", "bulk", "bulk_samples.qc.{}.expression_boxplots.png".format(quant_type)), dpi=300, bbox_inches="tight")

        # Pairwise correlations
        g = sns.clustermap(
            exp_matrix.corr(),
            row_cluster=True, col_cluster=True,
            xticklabels=True, yticklabels=True,
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        for item in g.ax_heatmap.get_xticklabels():
            item.set_rotation(90)
        g.fig.savefig(os.path.join("results", "bulk", "bulk_samples.qc.{}.correlation.png".format(quant_type)), dpi=300, bbox_inches="tight")

        # Heatmap and correlation on signature genes
        # derived from bulk
        # derived from scRNA
        for geneset in ["de_genes", "de_genes_bulk"]:
            g = sns.clustermap(
                exp_matrix.ix[eval(geneset)].dropna(),
                z_score=0,
                row_cluster=True, col_cluster=True,
                xticklabels=True, yticklabels=True,
                figsize=(15, 15))
            for item in g.ax_heatmap.get_yticklabels():
                item.set_rotation(0)
            for item in g.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
            g.fig.savefig(os.path.join("results", "bulk", "bulk_samples.qc.{}.{}.png".format(quant_type, geneset)), dpi=300, bbox_inches="tight")

            g = sns.clustermap(
                exp_matrix.ix[eval(geneset)].dropna().corr(),
                row_cluster=True, col_cluster=True,
                xticklabels=True, yticklabels=True,
                figsize=(15, 15))
            for item in g.ax_heatmap.get_yticklabels():
                item.set_rotation(0)
            for item in g.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
            g.fig.savefig(os.path.join("results", "bulk", "bulk_samples.qc.{}.{}.correlation.png".format(quant_type, geneset)), dpi=300, bbox_inches="tight")


def compare_bulk(matrix_norm, cond1="stimulated", cond2="unstimulated"):
    """
    Compare Bulk RNA-seq data with single-cell
    """
    from scipy.stats import pearsonr, spearmanr
    from sklearn.metrics import mean_squared_error
    from math import sqrt

    def normalize_quantiles_r(dataframe):
        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        rpy2.robjects.numpy2ri.activate()

        robjects.r('require("preprocessCore")')
        normq = robjects.r('normalize.quantiles')

        dataframe_norm = pd.DataFrame(
            np.array(normq(dataframe)),
            index=dataframe.index,
            columns=dataframe.columns
        )

        return dataframe_norm

    def seurat_regress(matrix_file, output_prefix):
        """
        """
        import rpy2.robjects as robj
        import pandas.rpy.common as com

        run = robj.r("""
            # module unload R
            # module load gcc/6.0.0
            # module load texlive
            # module load R/3.2.3
            # srun -J seuratBulk --mem 256000 -p develop -n 1 -c 8 --pty /cm/shared/apps/R/3.2.3/bin/R

            # matrix_file = "/home/arendeiro/projects/crop-seq/results/CROP-seq_Jurkat_TCR.count_matrix.gene_level.csv"
            # output_prefix = "/home/arendeiro/projects/crop-seq/results/CROP-seq_Jurkat_TCR.count_matrix.gene_level"

            # matrix_file = "/home/arendeiro/projects/crop-seq/results/CROP-seq_Jurkat_TCR.ESAT_count_matrix.csv"
            # output_prefix = "/home/arendeiro/projects/crop-seq/results/CROP-seq_Jurkat_TCR.ESAT_count_matrix"

            function(matrix_file, output_prefix){
                library(Seurat)
                library(dplyr)
                library(Matrix)

                counts <- read.csv(matrix_file, sep=",")

                # create experiment matrix
                colData = as.data.frame(t(counts[c(1, 2, 3), ]))
                colData$sample_name = rownames(colData)
                colnames(colData) = c("condition", "gene", "grna", "sample_name")
                colData = colData[-c(1), ]

                rownames(counts) = counts[, 1]
                countData = counts[-c(0, 1, 2, 3, 4), -c(1)]
                countData <- apply(countData, 2, function(x) {storage.mode(x) <- 'integer'; x})
                countData <- apply(countData, 2, function(x) {x[is.na(x)] <- 0; x})

                # Initialize the Seurat object with the raw data
                exp <- new("seurat", raw.data = countData)

                # Keep all genes expressed in at least 20 cells, keep all cells with >= 500 genes
                exp <- Setup(
                    exp, project="cropseq",
                    min.cells=5, min.genes=500,
                    do.logNormalize=T, total.expr=1e4)

                # Stash cell attributes
                exp <- AddMetaData(exp, colData$sample_name, "sample_name")
                exp <- AddMetaData(exp, colData$condition, "condition")
                exp <- AddMetaData(exp, colData$gene, "gene")
                exp <- AddMetaData(exp, colData$grna, "grna")

                # Calculate the percentage of mitochondrial genes and store it.
                mito.genes <- grep("^MT-", rownames(exp@data), value=T)
                percent.mito <- colSums(expm1(exp@data[mito.genes, ])) / colSums(expm1(exp@data))
                exp <- AddMetaData(exp, percent.mito, "percent.mito")
                # Calculate the percentage of ribosomal genes and store it.
                ribo.genes <- grep("^RP", rownames(exp@data), value=T)
                percent.ribo <- colSums(expm1(exp@data[ribo.genes, ])) / colSums(expm1(exp@data))
                exp <- AddMetaData(exp, percent.ribo, "percent.ribo")

                # Regress out
                exp_subset_regressed <- RegressOut(exp, latent.vars=c("nUMI", "percent.mito"))

                # Write as csv
                output = paste(output_prefix, "seurat_regressed.csv", sep=".")
                write.csv(exp_subset_regressed@scale.data, output, sep=",", quote=FALSE)
                return(as.matrix(exp_subset_regressed@data))
            }
        """)
        # convert to Python objects
        norm_matrix = com.convert_robj(
            run(matrix_file, output_prefix))

        return norm_matrix

    def rmse(x, y):
        return sqrt(mean_squared_error(x, y)), pd.np.nan

    def n_cells(x):
        return x.shape[1], pd.np.nan

    experiment = "CROP-seq_Jurkat_TCR"

    # get single-cell expression matrix and match grna names
    matrix_norm.columns = matrix_norm.columns.set_levels([matrix_norm.columns.levels[3].str.replace(".*library_", "")], level=[3])

    # Try using counts
    exp_assigned = pd.read_hdf(os.path.join(results_dir, "{}.digital_expression.{}genes.only_assigned.hdf5.gz".format(experiment, n_genes)), "exp_matrix", compression="gzip")
    exp_assigned = exp_assigned.T.reset_index()
    exp_assigned['replicate'] = exp_assigned['replicate'].astype(np.int64).astype(str)
    exp_assigned['gene'] = pd.np.nan
    exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 3, 'gene'] = pd.Series(exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 3, "grna"].str.split("_")).apply(lambda x: x[1])
    exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 4, 'gene'] = pd.Series(exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 4, "grna"].str.split("_")).apply(lambda x: x[2])
    exp_assigned.loc[exp_assigned['grna'].str.contains("CTRL"), 'gene'] = "CTRL"
    exp_assigned = exp_assigned.set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])
    exp_assigned = exp_assigned.T

    tpm = np.log2(1 + exp_assigned.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)
    tpm_quant = np.log2(1 + normalize_quantiles_r(exp_assigned.apply(lambda x: x / float(x.sum()), axis=0) * 1e4))

    # read gene expression matrix
    bitseq = pd.read_csv(os.path.join("results", "{}.count_matrix.gene_level.csv".format(experiment)), index_col=[0], header=range(4))
    esat = pd.read_csv(os.path.join("results", "{}.ESAT_count_matrix.csv".format(experiment)), index_col=[0], header=range(4)).fillna(0)

    bitseq_quant = np.log2(1 + normalize_quantiles_r(bitseq))
    esat_quant = np.log2(1 + normalize_quantiles_r(esat))

    bitseq = np.log2(1 + bitseq.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)
    esat = np.log2(1 + esat.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)

    # bitseq_seurat = seurat_regress(bitseq)
    # esat_seurat = seurat_regress(esat)
    bitseq_seurat = pd.read_csv(os.path.join("results", "{}.count_matrix.gene_level.seurat_regressed.csv".format(experiment)), index_col=[0])
    bitseq_seurat.columns = bitseq.columns
    esat_seurat = pd.read_csv(os.path.join("results", "{}.ESAT_count_matrix.seurat_regressed.csv".format(experiment)), index_col=[0])
    esat_seurat.columns = esat.columns

    # read in diff genes
    degs = pd.read_csv(
        os.path.join(results_dir, "{}.digital_expression.{}genes.scde.diff_expr.csv".format(experiment, n_genes)), index_col=0
    )
    degs = degs[~degs.index.str.contains("library|CTRL")]
    de_genes = degs[abs(degs["cZ"]) > 2].index.tolist()

    # Tasks:

    # 0). Compare expression in CTRL cells/samples only
    # 1). Compare expression levels across all genes
    # 2). Compare expression levels in signature genes
    # 2b). Compare expression in selected markers
    # 3). Compare signature (recall, agreement)
    # 4). Compare differential genes in KO vs CTRL (recall, agreement, enrichemnts)
    # x). Each single-cell vs Bulk
    single_quant_types = [
        ("seurat", matrix_norm), ("tpm", tpm), ("tpm_quant", tpm_quant)]

    bulk_quant_types = [
        ("bitseq", bitseq), ("esat", esat), ("bitseq_quant", bitseq_quant),
        ("esat_quant", esat_quant)]

    # 0).
    for i, condition in enumerate(matrix_norm.columns.get_level_values('condition').drop_duplicates()):
        fig, axis = plt.subplots(5, 5, figsize=(5 * 3, 5 * 3))
        axis = iter(axis.flatten())
        for gene_group, gene_filter in [("all_genes", tpm.index), ("de_genes", de_genes)]:
            for single_quant_type, single_exp_matrix in single_quant_types:
                for bulk_quant_type, bulk_exp_matrix in bulk_quant_types:
                    # remove gRNAs
                    single_exp_matrix = single_exp_matrix[~single_exp_matrix.index.str.contains("library|CTRL")]

                    # remove ribosomal, mitochondrial genes
                    single_exp_matrix = single_exp_matrix[(~single_exp_matrix.index.str.contains("^RP.*")) & (~single_exp_matrix.index.str.contains("^MT-"))]
                    bulk_exp_matrix = bulk_exp_matrix[(~bulk_exp_matrix.index.str.contains("^RP.*")) & (~bulk_exp_matrix.index.str.contains("^MT-"))]

                    # align indices
                    single_exp_matrix = single_exp_matrix.ix[bulk_exp_matrix.index].dropna()
                    bulk_exp_matrix = bulk_exp_matrix.ix[single_exp_matrix.index].dropna()
                    de_genes = [x for x in de_genes if x in bulk_exp_matrix.index]

                    a = single_exp_matrix[single_exp_matrix.columns[
                        (single_exp_matrix.columns.get_level_values('condition') == condition) &
                        (single_exp_matrix.columns.get_level_values('gene') == "CTRL")]].median(axis=1).ix[gene_filter].dropna()

                    b = bulk_exp_matrix[bulk_exp_matrix.columns[
                        (bulk_exp_matrix.columns.get_level_values('condition') == condition) &
                        (bulk_exp_matrix.columns.get_level_values('gene') == "CTRL")]].median(axis=1).ix[gene_filter].dropna()

                    ax = axis.next()
                    ax.scatter(a, b, alpha=0.75, s=3)
                    p = pearsonr(a, b)[0]
                    s = spearmanr(a, b)[0]
                    print(gene_group, single_quant_type, bulk_quant_type, p, s)
                    ax.text(a.min(), b.max(), "Pearson: {0:0.3f}\n Spearman: {1:0.3f}".format(p, s))
                    ax.set_title("{} {} {}".format(gene_group, single_quant_type, bulk_quant_type))
        sns.despine(fig)
        fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.scatter.CTRL_only.{}.png".format(condition)), dpi=300, bbox_inches="tight")

    # 1). & 2).
    # Select normalization method
    single_quant_type, single_exp_matrix = single_quant_types[2]
    bulk_quant_type, bulk_exp_matrix = bulk_quant_types[0]

    # remove gRNAs
    single_exp_matrix = single_exp_matrix[~single_exp_matrix.index.str.contains("library|CTRL")]

    # remove ribosomal, mitochondrial genes
    single_exp_matrix = single_exp_matrix[(~single_exp_matrix.index.str.contains("^RP.*")) & (~single_exp_matrix.index.str.contains("^MT-"))]
    bulk_exp_matrix = bulk_exp_matrix[(~bulk_exp_matrix.index.str.contains("^RP.*")) & (~bulk_exp_matrix.index.str.contains("^MT-"))]

    # align indices
    single_exp_matrix = single_exp_matrix.ix[bulk_exp_matrix.index].dropna()
    bulk_exp_matrix = bulk_exp_matrix.ix[single_exp_matrix.index].dropna()
    de_genes = [x for x in de_genes if x in bulk_exp_matrix.index]

    # get single-cell expression matrix and match grna names
    single_exp_matrix.columns = single_exp_matrix.columns.set_levels([single_exp_matrix.columns.levels[3].str.replace(".*library_", "")], level=[3])

    # Run
    comparisons = pd.DataFrame()
    for gene_group, gene_filter in [("all_genes", bulk_exp_matrix.index), ("de_genes", de_genes)]:
        # At gene level or at grna level
        for level in ["gene", "grna"]:
            fig, axis = plt.subplots(
                len(single_exp_matrix.columns.get_level_values('condition').drop_duplicates()),
                len(single_exp_matrix.columns.get_level_values(level).drop_duplicates()),
                figsize=(
                    4 * len(single_exp_matrix.columns.get_level_values(level).drop_duplicates()),
                    4 * len(single_exp_matrix.columns.get_level_values('condition').drop_duplicates())
                ), sharex=True, sharey=True
            )
            q = pd.DataFrame()
            for i, condition in enumerate(single_exp_matrix.columns.get_level_values('condition').drop_duplicates()):
                # Compare within each knockout
                for j, ko in enumerate(single_exp_matrix.columns.get_level_values(level).drop_duplicates()):
                    print(bulk_quant_type, gene_group, level, condition, ko)

                    s = single_exp_matrix.columns[(single_exp_matrix.columns.get_level_values("condition") == condition) & (single_exp_matrix.columns.get_level_values(level) == ko)]
                    b = bulk_exp_matrix.columns[(bulk_exp_matrix.columns.get_level_values("condition") == condition) & (bulk_exp_matrix.columns.get_level_values(level) == ko)]

                    if b.shape[0] == 0 or s.shape[0] == 0:
                        continue

                    single_cell = single_exp_matrix[s]
                    sm = single_cell.median(axis=1).ix[gene_filter]
                    bulk = bulk_exp_matrix[b]
                    bm = bulk.median(axis=1).ix[gene_filter]

                    q['{} {} {}'.format("single_cell", condition, ko)] = sm
                    q['{} {} {}'.format("bulk", condition, ko)] = bm

                    # Mean of cells vs Bulk
                    res = list()
                    for metric in [pearsonr, spearmanr, rmse, n_cells]:
                        if metric != n_cells:
                            r = metric(sm, bm)
                        else:
                            r = metric(single_cell)
                        res.append(r)

                        comparisons = comparisons.append(pd.Series(
                            [gene_group, condition, level, ko, metric.__name__] + list(r),
                            index=["gene_group", "condition", "level", "KO", "metric", "value", "p_value"]), ignore_index=True)

                    axis[i][j].scatter(bm, sm, alpha=0.75, s=4)

                    axis[i][j].set_title("{} {}".format(condition, ko))
                    axis[i][j].text(bm.max(), sm.max(), "Pearson: {0:0.3f}\n Spearman: {1:0.3f}\nRMSE: {2:0.3f}".format(*[x[0] for x in res]))
                    axis[i][j].set_xlabel("Bulk 3' RNA-seq")
                    axis[i][j].set_ylabel("Mean of single-cells")
            sns.despine(fig)
            fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.scatter.{}.png".format(level)), dpi=300, bbox_inches="tight")

            comparisons.to_csv(os.path.join("results", "bulk", "bulk_single-cell_comparison.metrics.csv"))

            comparisons2 = comparisons.copy()
            comparisons2.loc[comparisons2['metric'] == "n_cells", "value"] = np.log10(comparisons2.loc[comparisons2['metric'] == "n_cells", "value"])

            # heatmap with metrics
            fig, axis = plt.subplots(
                len(comparisons2["metric"].drop_duplicates()),
                len(single_exp_matrix.columns.get_level_values('condition').drop_duplicates()),
            )
            for i, condition in enumerate(single_exp_matrix.columns.get_level_values('condition').drop_duplicates()):
                for j, metric in enumerate(comparisons2["metric"].drop_duplicates()):
                    pivot = pd.pivot_table(comparisons2[
                        (comparisons2["gene_group"] == gene_group) &
                        (comparisons2["level"] == level) &
                        (comparisons2["condition"] == condition) &
                        (comparisons2["metric"] == metric)], values="value", index="metric", columns="KO")
                    sns.heatmap(pivot, ax=axis[j][i])

                    if j != 3:
                        axis[j][i].set_xticklabels(axis[j][i].get_xticklabels(), visible=False)
                        axis[j][i].set_xlabel(None, visible=False)
                    else:
                        axis[j][i].set_xticklabels(axis[j][i].get_xticklabels(), rotation=90)
            fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.{}.metric_heatmap.{}.png".format(gene_group, level)), dpi=300, bbox_inches="tight")

            # scatter of n_cells vs metrics
            pivot = pd.pivot_table(comparisons2[
                (comparisons2["gene_group"] == gene_group) &
                (comparisons2["level"] == level)], values="value", index=["condition", "KO"], columns="metric")
            g = sns.pairplot(data=pivot.reset_index().set_index("KO"), hue="condition")
            g.fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.{}.metric_scatter.{}.png".format(gene_group, level)), dpi=300, bbox_inches="tight")

            comparisons.to_csv(os.path.join("results", "bulk", "bulk_single-cell_comparison.metrics.csv"))

    # 2b).
    # CD69, CD82, etc...
    # At gene level or at grna level
    for level in ["gene", "grna"]:
        print(level)
        markers = ["CD69", "CD82", "PDCD1", "CD38", "BCL7A", "CDC20", "TUBB", "ADA", "TUBA1B"]

        sm = single_exp_matrix.T.groupby(level=['condition', level]).mean()[markers]
        sm["data_type"] = "single_cell"
        bm = bulk_exp_matrix.T.groupby(level=['condition', level]).mean()[markers]
        # bm = bm.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
        bm["data_type"] = "bulk"
        bm = bm.ix[sm.index].dropna()
        sm = sm.ix[bm.index].dropna()
        e = pd.melt(sm.reset_index(), id_vars=['condition', level, 'data_type']).append(pd.melt(bm.reset_index(), id_vars=['condition', level, 'data_type']))

        g = sns.FacetGrid(data=e, col="gene_name", row="data_type", hue="condition", sharey=False, sharex=True)
        g.map(sns.stripplot, level, "value")
        g.add_legend()
        for ax in g.axes.flat:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.marker_genes.stripplot.{}.png".format(level)), dpi=300, bbox_inches="tight")

        g = sns.FacetGrid(data=e, col="gene_name", row="condition", hue="data_type", sharey=False, sharex=True)
        g.map(sns.stripplot, level, "value")
        g.add_legend()
        for ax in g.axes.flat:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.marker_genes.stripplot.techniq_together.{}.png".format(level)), dpi=300, bbox_inches="tight")

        fig, axis = plt.subplots(
            len(single_exp_matrix.columns.get_level_values('condition').drop_duplicates()),
            len(markers),
            figsize=(
                4 * len(markers),
                4 * len(single_exp_matrix.columns.get_level_values('condition').drop_duplicates()))
        )
        for i, condition in enumerate(single_exp_matrix.columns.get_level_values('condition').drop_duplicates()):
            for j, marker in enumerate(markers):
                axis[i][j].scatter(sm.loc[condition, marker], bm.loc[condition, marker])
                axis[i][j].set_title("{} {}".format(condition, marker))
                c, p = pearsonr(sm.loc[condition, marker], bm.loc[condition, marker])
                axis[i][j].text(max(sm.loc[condition, marker]), max(bm.loc[condition, marker]), "r = {0:0.3f}\np = {1:0.3f}".format(c, p))
        sns.despine(fig)
        fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.marker_genes.scatter.{}.png".format(level)), dpi=300, bbox_inches="tight")

    # 3).

    # 4).
    def agreement(x, y):
        overlap = x.index.isin(y.index).sum()
        recall = overlap / float(y.shape[0])
        return overlap, recall

    performance = pd.DataFrame()
    for i, condition in enumerate(bitseq.columns.levels[1]):
        for j, ko in enumerate(bitseq.columns.levels[2]):
            print(condition, ko)
            # get differential
            # read in single
            try:
                single_diff = pd.read_csv(os.path.join("results", "CROP-seq_Jurkat_TCR.digital_expression.500genes.scde.diff_expr.{}.{}.csv".format(condition, ko)))
            except IOError:
                continue
            single_diff = single_diff[~single_diff.index.str.contains("library|CTRL")]
            single_diff["abs_cZ"] = abs(single_diff["cZ"])

            # read in Bulk
            try:
                bulk_diff = pd.read_csv(os.path.join("results", "bulk", "deseq_bitseq_{}_{}.diff.csv".format(condition, ko)))
            except IOError:
                continue

            # remove gRNAs
            single_diff = single_diff[~single_diff.index.str.contains("library|CTRL")]

            # remove ribosomal, mitochondrial genes
            single_diff = single_diff[(~single_diff.index.str.contains("^RP.*")) & (~single_diff.index.str.contains("^MT-"))]
            bulk_diff = bulk_diff[(~bulk_diff.index.str.contains("^RP.*")) & (~bulk_diff.index.str.contains("^MT-"))]

            # match the indices
            single_diff = single_diff.ix[bulk_diff.index].dropna()
            bulk_diff = bulk_diff.ix[single_diff.index].dropna()

            # sort by best
            single_diff = single_diff.sort_values("abs_cZ", ascending=False)
            bulk_diff = bulk_diff.sort_values("padj", ascending=True)

            # compute recall (based on fixed thresholds)
            for method in ["top250", "pvalue"]:
                x = single_diff[abs(single_diff['cZ']) > 2]
                if method == "top250":
                    gold_y = bulk_diff.head(250)
                elif method == "pvalue":
                    gold_y = bulk_diff[bulk_diff['padj'] < 0.1]
                r, o = agreement(x, gold_y)
                performance = performance.append(pd.Series(
                    [condition, ko, "std", method, x.shape[0], gold_y.shape[0], o, r],
                    index=["condition", "ko", "threshold", "method", "sc_size", "bulk_size", "overlap", "recall"]), ignore_index=True)

                # compute agreement (on a ranked gradient)
                for g in np.arange(0.01, 1.01, 0.01):

                    x = single_diff.irow(range(0, int(np.round(bulk_diff.shape[0] * g))))
                    y = bulk_diff.irow(range(0, int(np.round(bulk_diff.shape[0] * g))))

                    o, r = agreement(x, y)
                    oo, rr = agreement(x, gold_y)
                    performance = performance.append(pd.Series(
                        [condition, ko, g, method, x.shape[0], y.shape[0], o, r, oo, rr],
                        index=["condition", "ko", "threshold", "method", "sc_size", "bulk_size", "overlap", "recall", "gold_overlap", "gold_recall"]), ignore_index=True)
    performance.to_csv(os.path.join("results", "bulk", "bulk_single-cell_comparison.diff_genes.performance.csv"))

    for i, condition in enumerate(bitseq.columns.levels[1]):
        for method in ["top250", "pvalue"]:
            # Plot recall of both
            pivot = pd.pivot_table(performance[
                (performance["condition"] == condition) & (performance["threshold"] != "std") & (performance["method"] == method)
            ], index="ko", columns="threshold", values="recall")
            fig, axis = plt.subplots(1)
            sns.heatmap(pivot, ax=axis)
            axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
            axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
            fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.diff_genes.recall_per_rank_percentile.{}.{}.png".format(condition, method)), dpi=300, bbox_inches="tight")

            # Plot recall of gold standard
            pivot = pd.pivot_table(performance[
                (performance["condition"] == condition) & (performance["threshold"] != "std") & (performance["method"] == method)], index="ko", columns="threshold", values="gold_recall")
            fig, axis = plt.subplots(1)
            sns.heatmap(pivot, ax=axis, vmin=0, vmax=1)
            axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
            axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
            fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.diff_genes.recall_of_gold_per_rank_percentile.{}.{}.png".format(condition, method)), dpi=300, bbox_inches="tight")

    # Get enrichments
    n_top_genes = 500
    # save text file with gene names
    for i, condition in enumerate(bitseq.columns.levels[1]):
        for j, ko in enumerate(bitseq.columns.levels[2]):
            print(condition, ko)
            # get differential
            # read in single
            try:
                single_diff = pd.read_csv(os.path.join("results", "CROP-seq_Jurkat_TCR.digital_expression.500genes.scde.diff_expr.{}.{}.csv".format(condition, ko)))
            except IOError:
                continue

            # read in Bulk
            try:
                bulk_diff = pd.read_csv(os.path.join("results", "bulk", "deseq_bitseq_{}_{}.diff.csv".format(condition, ko)))
            except IOError:
                continue

            # remove gRNAs
            single_diff = single_diff[~single_diff.index.str.contains("library|CTRL")]

            # remove ribosomal, mitochondrial genes
            single_diff = single_diff[(~single_diff.index.str.contains("^RP.*")) & (~single_diff.index.str.contains("^MT-"))]
            bulk_diff = bulk_diff[(~bulk_diff.index.str.contains("^RP.*")) & (~bulk_diff.index.str.contains("^MT-"))]

            # sort by best
            single_diff["abs_cZ"] = abs(single_diff["cZ"])
            single_diff = single_diff.sort_values("abs_cZ", ascending=False)
            bulk_diff = bulk_diff.sort_values("padj", ascending=True)

            # get top 250
            single_diff.head(n_top_genes).index.to_series().to_csv(os.path.join("results", "bulk", "enrichr", "single_cell.{}.{}.gene_symbols.txt".format(condition, ko)), index=False)
            bulk_diff.head(n_top_genes).index.to_series().to_csv(os.path.join("results", "bulk", "enrichr", "bulk.{}.{}.gene_symbols.txt".format(condition, ko)), index=False)

            # submit to Enrichr
            for t in ["bulk", "single_cell"]:
                cmd = "sbatch -J {} -o {} -p shortq --mem 8000 -c 1 ~/jobs/run_Enrichr.sh {} {}".format(
                    "enrichr.{}.{}".format(condition, ko), os.path.join("results", "bulk", "enrichr", "{}.{}.{}.log".format(t, condition, ko)),
                    os.path.join("results", "bulk", "enrichr", "{}.{}.{}.gene_symbols.txt".format(t, condition, ko)),
                    os.path.join("results", "bulk", "enrichr", "{}.{}.{}.enrichr.csv".format(t, condition, ko))
                )
                os.system(cmd)

    # collect
    enrichments = pd.DataFrame()
    for i, condition in enumerate(bitseq.columns.levels[1]):
        for j, ko in enumerate(bitseq.columns.levels[2]):
            print(condition, ko)
            # get differential
            # read in single
            try:
                single_diff = pd.read_csv(os.path.join("results", "bulk", "enrichr", "single_cell.{}.{}.enrichr.csv".format(condition, ko)))
            except IOError:
                continue

            # read in Bulk
            try:
                bulk_diff = pd.read_csv(os.path.join("results", "bulk", "enrichr", "bulk.{}.{}.enrichr.csv".format(condition, ko)))
            except IOError:
                continue

            single_diff["condition"] = bulk_diff["condition"] = condition
            single_diff["ko"] = bulk_diff["ko"] = ko
            single_diff["data_type"] = "single_cell"
            bulk_diff["data_type"] = "bulk"

            enrichments = enrichments.append(single_diff).append(bulk_diff)

    enrichments.to_csv(os.path.join("results", "bulk", "enrichr", "all_comparisons.{}topgenes.enrichr.csv".format(n_top_genes)), index=False)

    for i, condition in enumerate(bitseq.columns.levels[1]):
        for gene_set_library in enrichments['gene_set_library'].drop_duplicates():
            print(condition, gene_set_library)
            # pivot
            pivot = -np.log10(pd.pivot_table(
                enrichments[(enrichments['condition'] == condition) & (enrichments['gene_set_library'] == gene_set_library)],
                index=['data_type', 'ko'], columns="description", values="p_value").fillna(1))

            variable = (pivot.std() / pivot.sum()).sort_values()
            try:
                # plot
                g = sns.clustermap(
                    pivot[variable.tail(50).index],
                    row_cluster=True, col_cluster=True,
                    xticklabels=True, yticklabels=True,
                    figsize=(15, 15))
                for item in g.ax_heatmap.get_yticklabels():
                    item.set_rotation(0)
                for item in g.ax_heatmap.get_xticklabels():
                    item.set_rotation(90)
                g.savefig(os.path.join("results", "bulk", "enrichr", "enrichr.{}.{}.clustermap.png".format(condition, gene_set_library)), dpi=300, bbox_inches="tight")
                g = sns.clustermap(
                    pivot[variable.tail(50).index].T.corr(),
                    row_cluster=True, col_cluster=True,
                    xticklabels=True, yticklabels=True,
                    figsize=(15, 15))
                for item in g.ax_heatmap.get_yticklabels():
                    item.set_rotation(0)
                for item in g.ax_heatmap.get_xticklabels():
                    item.set_rotation(90)
                g.savefig(os.path.join("results", "bulk", "enrichr", "enrichr.{}.{}.clustermap.correlation.png".format(condition, gene_set_library)), dpi=300, bbox_inches="tight")
            except FloatingPointError:
                continue


def flow_analysis():
    """
    """
    from FlowCytometryTools import FCMeasurement
    from FlowCytometryTools import ThresholdGate, PolyGate
    from scipy.stats import pearsonr

    flow = pd.read_csv(os.path.join("metadata", "flow_analysis.csv"))

    # parse fcs files,
    # gate on single, live cells
    # extract cells passing crteria
    fig0, axis0 = plt.subplots(8, 12, sharex=True, sharey=True, figsize=(12 * 2, 8 * 2))
    fig1, axis1 = plt.subplots(8, 12, sharex=True, sharey=True, figsize=(12 * 2, 8 * 2))
    fig2, axis2 = plt.subplots(8, 12, sharex=True, sharey=True, figsize=(12 * 2, 8 * 2))
    fcs_cells = dict()
    for i, fcs_file in enumerate(flow['fcs_file']):
        if (flow.loc[flow['fcs_file'] == fcs_file, "failed"]).squeeze() is True:
            continue
        sample = FCMeasurement(ID=fcs_file, datafile=os.path.join("flow_data", fcs_file))

        # extract measurements:
        # 1st gating
        gate1 = ThresholdGate(17500, 'FSC-A', region='above')
        gate2 = ThresholdGate(100000, 'SSC-A', region='below')
        gate3 = ThresholdGate(200000, 'FSC-A', region='below')
        gate4 = ThresholdGate(7000, 'SSC-A', region='above')
        _ = sample.plot(['FSC-A', 'SSC-A'], bins=200, ax=axis0.flat[i], gates=[gate1, gate2, gate3, gate4])

        for gate in [gate1, gate2, gate3, gate4]:
            sample = sample.gate(gate, apply_now=True)

        # to plot after gating
        # _ = sample.plot(['FSC-A', 'SSC-A'], bins=2000, gates=[gate1, gate2, gate3, gate4])

        # 2nd gate
        gate5 = PolyGate([
            (1.782e+04, 1.896e+04), (1.060e+05, 1.022e+05),
            (1.898e+05, 1.611e+05), (1.990e+05, 1.609e+05), (1.994e+05, 1.168e+05),
            (1.060e+05, 5.667e+04), (4.525e+04, 2.380e+04), (1.801e+04, 1.129e+04),
            (1.782e+04, 1.875e+04), (1.856e+04, 1.916e+04)], ('FSC-A', 'FSC-H'), region='in', name='FSC-A_FSC-H')
        _ = sample.plot(['FSC-A', 'FSC-H'], bins=200, ax=axis1.flat[i], gates=[gate5])

        sample = sample.gate(gate5, apply_now=True)

        # to plot after gating
        # _ = sample.plot(['FSC-A', 'FSC-H'], bins=2000, gates=[gate5])

        # 3rd gate
        gate6 = ThresholdGate(75, 'R/A APC-Cy7-A', region='below')
        _ = sample.plot(['FSC-A', 'R/A APC-Cy7-A'], bins=200, ax=axis2.flat[i], gates=[gate6])
        sample = sample.gate(gate6, apply_now=True)

        # extract good cells
        data = sample.data

        fcs_cells[fcs_file] = data.drop("Time", axis=1)

        axis0.flat[i].set_xlim((0, 250e3))
        axis0.flat[i].set_xlim((0, 150000))
        condition = (flow.loc[flow['fcs_file'] == fcs_file, "condition"]).squeeze()
        grna = (flow.loc[flow['fcs_file'] == fcs_file, "grna"]).squeeze()
        axis0.flat[i].set_title("{} - {}".format(condition, grna))
        axis1.flat[i].set_title("{} - {}".format(condition, grna))
        axis2.flat[i].set_title("{} - {}".format(condition, grna))

    for axis in [axis0, axis1, axis2]:
        for a in axis.flat:
            a.set_xlabel(None, visible=False)
            a.set_ylabel(None, visible=False)
            a.set_xticklabels(a.get_xticklabels(), visible=False)
            a.set_yticklabels(a.get_yticklabels(), visible=False)
    for f in [fig0, fig1, fig2]:
        sns.despine(f)
    fig0.savefig(os.path.join("results", "flow", "flow.all_samples.FSCA_SSCA.png"), bbox_inches="tight", dpi=300)
    fig1.savefig(os.path.join("results", "flow", "flow.all_samples.FSCA_FSCH.png"), bbox_inches="tight", dpi=300)
    fig2.savefig(os.path.join("results", "flow", "flow.all_samples.FSCA_LD.png"), bbox_inches="tight", dpi=300)

    #

    # Plot pregated mean of positive cells vs marker gene expression per gRNA
    marker_gene_mapping = {
        "CD38": "CD38", "CD69": "CD69", "PD1": "PDCD1", "pRelA": "RELA", "CD154": "CD40LG", "CD25": "IL2RA"
    }
    matrix_norm_t = matrix_norm.T

    fig, axis = plt.subplots(
        len(flow['condition'].drop_duplicates()), len(marker_gene_mapping.keys()), sharex=False, sharey=False,
        figsize=(
            len(marker_gene_mapping.keys()) * 4,
            len(flow['condition'].drop_duplicates()) * 4))

    gene_colors = dict(zip(flow['gene'].drop_duplicates(), sns.color_palette("colorblind") * 10))
    cm = plt.get_cmap('gist_rainbow')
    gene_colors = dict(zip(flow['gene'].drop_duplicates(), [cm(1. * i / len(flow['gene'].drop_duplicates())) for i in range(len(flow['gene'].drop_duplicates()))]))
    for i, condition in enumerate(flow['condition'].drop_duplicates()):
        for j, marker in enumerate(marker_gene_mapping.keys()):
            facs_means = dict()
            scRNA_means = dict()
            for grna in flow[flow['condition'] == condition]['grna'].drop_duplicates():
                gene = (flow.loc[(flow['condition'] == condition) & (flow['grna'] == grna), "gene"]).squeeze()
                facs_mean = flow[(flow['condition'] == condition) & (flow['grna'] == grna)]["%{}+".format(marker)].squeeze()
                scRNA_mean = matrix_norm_t.ix[
                    matrix_norm_t.index[
                        (matrix_norm_t.index.get_level_values("condition") == condition) &
                        (matrix_norm_t.index.get_level_values("grna") == ("Tcrlibrary_" + grna if not grna.startswith("CTRL") else grna))
                    ]][marker_gene_mapping[marker]].mean()

                axis[i][j].scatter(scRNA_mean, facs_mean, s=10, alpha=0.7, color=gene_colors[gene])
                axis[i][j].set_title("{} {}".format(condition, marker))

                facs_means[grna] = facs_mean
                scRNA_means[grna] = scRNA_mean

            c, p = pearsonr(scRNA_means.values(), facs_means.values())
            axis[i][j].text(max(scRNA_means.values()), max(facs_means.values()), "r = {0:0.3f}\np = {1:0.3f}".format(c, p))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "flow", "flow.all_samples.mean_flow_mean_scRNA.grna.png"), bbox_inches="tight", dpi=300)

    # Only CD69
    fig, axis = plt.subplots(
        len(flow['condition'].drop_duplicates()), 1, sharex=False, sharey=False,
        figsize=(
            1 * 4,
            len(flow['condition'].drop_duplicates()) * 4))

    gene_colors = dict(zip(flow['gene'].drop_duplicates(), sns.color_palette("colorblind") * 10))
    cm = plt.get_cmap('gist_rainbow')
    gene_colors = dict(zip(flow['gene'].drop_duplicates(), [cm(1. * i / len(flow['gene'].drop_duplicates())) for i in range(len(flow['gene'].drop_duplicates()))]))
    for i, condition in enumerate(flow['condition'].drop_duplicates()):
        for j, marker in enumerate(marker_gene_mapping.keys()):
            if marker != "CD69":
                continue
            facs_means = dict()
            scRNA_means = dict()
            for grna in flow[flow['condition'] == condition]['grna'].drop_duplicates():
                gene = (flow.loc[(flow['condition'] == condition) & (flow['grna'] == grna), "gene"]).squeeze()
                facs_mean = flow[(flow['condition'] == condition) & (flow['grna'] == grna)]["%{}+".format(marker)].squeeze()
                scRNA_mean = matrix_norm_t.ix[
                    matrix_norm_t.index[
                        (matrix_norm_t.index.get_level_values("condition") == condition) &
                        (matrix_norm_t.index.get_level_values("grna") == ("Tcrlibrary_" + grna if not grna.startswith("CTRL") else grna))
                    ]][marker_gene_mapping[marker]].mean()

                axis[i].scatter(scRNA_mean, facs_mean, s=10, alpha=0.7, color=gene_colors[gene])
                axis[i].text(scRNA_mean, facs_mean, grna)
                axis[i].set_title("{} {}".format(condition, marker))

                facs_means[grna] = facs_mean
                scRNA_means[grna] = scRNA_mean

            c, p = pearsonr(scRNA_means.values(), facs_means.values())
            axis[i].text(max(scRNA_means.values()), max(facs_means.values()), "r = {0:0.3f}\np = {1:0.3f}".format(c, p))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "flow", "flow.all_samples.mean_flow_mean_scRNA.grna.CD69_only.png"), bbox_inches="tight", dpi=300)

    # Plot pregated mean of positive cells vs marker gene expression per gene
    fig, axis = plt.subplots(
        len(flow['condition'].drop_duplicates()), len(marker_gene_mapping.keys()), sharex=False, sharey=False,
        figsize=(
            len(marker_gene_mapping.keys()) * 4,
            len(flow['condition'].drop_duplicates()) * 4))

    gene_colors = dict(zip(flow['gene'].drop_duplicates(), sns.color_palette("colorblind") * 10))
    cm = plt.get_cmap('gist_rainbow')
    gene_colors = dict(zip(flow['gene'].drop_duplicates(), [cm(1. * i / len(flow['gene'].drop_duplicates())) for i in range(len(flow['gene'].drop_duplicates()))]))
    for i, condition in enumerate(flow['condition'].drop_duplicates()):
        for j, marker in enumerate(marker_gene_mapping.keys()):
            facs_means = dict()
            scRNA_means = dict()
            for gene in flow[flow['condition'] == condition]['gene'].drop_duplicates():
                facs_mean = flow[(flow['condition'] == condition) & (flow['gene'] == gene)]["%{}+".format(marker)].squeeze().mean()
                scRNA_mean = matrix_norm_t.ix[
                    matrix_norm_t.index[
                        (matrix_norm_t.index.get_level_values("condition") == condition) &
                        (matrix_norm_t.index.get_level_values("gene") == gene)
                    ]][marker_gene_mapping[marker]].mean()

                axis[i][j].scatter(scRNA_mean, facs_mean, s=10, alpha=0.7, color=gene_colors[gene])
                axis[i][j].set_title("{} {}".format(condition, marker))

                facs_means[gene] = facs_mean
                scRNA_means[gene] = scRNA_mean

            c, p = pearsonr(scRNA_means.values(), facs_means.values())
            axis[i][j].text(max(scRNA_means.values()), max(facs_means.values()), "r = {0:0.3f}\np = {1:0.3f}".format(c, p))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "flow", "flow.all_samples.mean_flow_mean_scRNA.gene.png"), bbox_inches="tight", dpi=300)


prj = Project(os.path.join("metadata", "config.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = results_dir = os.path.join("results")

sample_annotation = prj.sheet.df

# get guide annotation
guide_annotation = os.path.join("metadata", "guide_annotation.csv")
guide_annotation = pd.read_csv(guide_annotation)

n_genes = 500
experiment = prj.sheet.df['experiment'].drop_duplicates()[4]


# get expression
for n_genes in [500]:
    for experiment in prj.sheet.df['experiment'].drop_duplicates():
        print(experiment)

        counts_file = os.path.join(results_dir, "{}.digital_expression.{}genes.only_assigned.hdf5.gz".format(experiment, n_genes))
        exp_assigned = pd.read_hdf(counts_file, "exp_matrix", compression="gzip")
        exp_assigned = exp_assigned.T.reset_index()
        exp_assigned['replicate'] = exp_assigned['replicate'].astype(np.int64).astype(str)
        exp_assigned['gene'] = pd.np.nan
        exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 3, 'gene'] = pd.Series(exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 3, "grna"].str.split("_")).apply(lambda x: x[1])
        exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 4, 'gene'] = pd.Series(exp_assigned.loc[pd.Series(exp_assigned['grna'].str.split("_")).apply(len) == 4, "grna"].str.split("_")).apply(lambda x: x[2])
        exp_assigned.loc[exp_assigned['grna'].str.contains("CTRL"), 'gene'] = "CTRL"
        exp_assigned = exp_assigned.set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])
        exp_assigned = exp_assigned.T

        # exp_assigned.columns = exp_assigned.columns.set_levels([exp_assigned.columns.levels[:-1], exp_assigned.columns.levels[-1].astype(str)])

        # Normalize
        # Approach 1:
        # normalize by total
        matrix_norm = normalize(exp_assigned, experiment=experiment, kind="total")

        # Approach 2:
        # regress out based on total number and MT-genes using Seurat
        # matrix_norm = normalize(counts_file, kind="seurat")

        # hdf5_file = os.path.join(results_dir, "{}.digital_expression.{}genes.only_assigned.seurat_regressed.hdf5.gz".format(experiment, n_genes))
        # matrix_norm = read_seurat_hdf5(hdf5_file)

        # Filtering/Cleanup
        # remove lowly expressed genes
        df = matrix_norm[matrix_norm.sum(1) > 10]

        # remove gRNAs
        df = df[~df.index.str.contains("library|CTRL")]

        # remove ribosomal, mitochondrial genes
        df = df[~((df.index.str.contains("^RP")) | (df.index.str.contains("^MRP")) | (df.index.str.contains("^MT-")))]

        # remove Essential genes
        df = df[df.columns[~df.columns.get_level_values("gene").isin(["DHODH", "MVD", "TUBB"])]]

        # Unsupervised

        # Approach 1:
        # apply dimentionality reduction methods/clustering
        # and discover biological axis related with stimulation
        unsupervised(df)
        differential_genes(df, experiment=experiment, method="pca")

        # Approach 2 (just a bit supervised though):
        # get differential genes between conditions from CTRL cells only
        # use this signature to position each cell
        # observe deviation of groups of cells with gRNA targeting the same
        # gene coming from either condition to the mean

        # variant A:

        # get differential genes with some method
        # differential_genes(df, experiment=experiment, method="scde")  # <- todo: add scde R code wrapped here

        # visualize signature and make signature position assignments per cell/gRNA/gene
        stimulation_signature(matrix_norm, method="pca")

        # Compare with FACS measurements
        inspect_bulk()
        compare_bulk(matrix_norm)
        flow_analysis()

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
        gather_scde()
        # explore_knockouts_scde()

        # Plot the difference between stimulated/unstimulated for same genes

        # B.
        # Get enrichments for all DEGs from the one vs all or one vs CTRL comparisons.
        # Plot knockout vs enrichment
