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
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rc('text', usetex=False)

import sys
sys.setrecursionlimit(10000)


def normalize(df, experiment="", kind="total", results_dir="results"):
    def normalize_by_total(df):
        """
        Normalize expression by total number of counts per cell.
        """
        return np.log2(1 + df.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)

    if kind == "total":
        norm = normalize_by_total(df)
        norm.to_hdf(os.path.join(results_dir, "digital_expression.500genes.{}.log2_tpm.hdf5.gz".format(experiment)), "log2_tpm", compression="gzip")
        return norm


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


def unsupervised(df, experiment="CROP-seq_Jurkat_TCR", filter_low=False, results_dir="results"):
    # Inspect
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE, MDS, LocallyLinearEmbedding, SpectralEmbedding, Isomap

    methods = ["PCA", "TSNE", "LocallyLinearEmbedding", "SpectralEmbedding", "Isomap", "MDS"]

    df2 = df[df.sum().sort_values().tail(int(df.shape[1] * 0.25)).index]  # get top 25% covered cells

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

    # Expression of some markers in PCA
    pca = PCA(n_components=2)
    fitted = pca.fit(df2.T).transform(df2.T)
    # eigenvalues = pca.explained_variance_ratio_
    # loadings = pca.components_
    # z_fitted = np.dot(fitted, loadings)

    fitted = pd.DataFrame(fitted, index=df2.columns, columns=["PC1", "PC2"])

    fig, axis = plt.subplots(1, 5, figsize=(4 * 5, 4 * 1))
    for level in range(2):
        axis[level].scatter(fitted["PC1"], fitted["PC2"], s=10, alpha=0.4, color=get_level_colors(df2.columns)[level])
    cmap = matplotlib.cm.Blues
    axis[5 - 1].scatter(fitted["PC1"], fitted["PC2"], s=15, alpha=0.5, color=cmap(z_score(df2.ix["TRAC"])))
    axis[5 - 2].scatter(fitted["PC1"], fitted["PC2"], s=15, alpha=0.5, color=cmap(z_score(df2.ix["JARID2"])))
    axis[5 - 3].scatter(fitted["PC1"], fitted["PC2"], s=15, alpha=0.5, color=cmap(z_score(df2.ix["PCNA"])))
    fig.savefig(os.path.join(results_dir, "{}.PCA_clustering.specific_genes.svg".format(experiment)), bbox_inches="tight")


def significant_perturbation(df, df_bulk, de_genes, results_dir="results", experiment="CROP-seq_Jurkat_TCR"):
    """
    Assess whether a gRNA perturbation is significant.
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.sandbox.stats.multicomp import multipletests
    from statsmodels.nonparametric.smoothers_lowess import lowess
    from scipy.stats import norm
    from scipy.stats import combine_pvalues

    def z_score(x):
        return (x - x.mean()) / x.std()

    # Calculate perturbation p-values
    results = pd.DataFrame()
    for data_type, tmp_df in [("CROP", df), ("Bulk", df_bulk)]:
        for subset, index in [("all_genes", df.index), ("sig_genes", de_genes)]:

            # groupby gRNA, reduce to median
            df_guide = tmp_df.ix[index].dropna().T.groupby(level=["condition", "gene", "grna"]).median().T

            for condition in df_guide.columns.levels[0]:
                ctrl = df_guide.loc[
                    :, (df_guide.columns.get_level_values('grna').str.contains("CTRL")) & (df_guide.columns.get_level_values('condition') == condition)
                ].median(axis=1)

                # fit gaussian on control gRNA expression (null)
                params = norm.fit(z_score(ctrl))

                for grna in df_guide.columns.levels[2]:
                    print(data_type, subset, condition, grna)
                    g = df_guide.loc[
                        :, (df_guide.columns.get_level_values('grna') == grna) & (df_guide.columns.get_level_values('condition') == condition)
                    ].squeeze()
                    n_cells = tmp_df.columns[
                        (tmp_df.columns.get_level_values('grna') == grna) & (tmp_df.columns.get_level_values('condition') == condition)
                    ].shape[0]

                    if g.empty or ctrl.empty:
                        continue

                    # get MannWhitney p-value
                    s, p = mannwhitneyu(g.astype(np.float128), ctrl.astype(np.float128))

                    # get p-value from fitted dist
                    es, ep = combine_pvalues(norm.sf(abs(z_score(g)), *params) * 2)  # twosided p-value

                    results = results.append(pd.Series(
                        [data_type, subset, condition, g.name[1], grna, s, p, n_cells, es, ep],
                        index=["data_type", "subset", "condition", "gene", "grna", "stat", "p_value", "n_cells", "estat", "ep_value"]
                    ), ignore_index=True)

    results = results[~results["grna"].str.contains("CTRL")]

    # correct p-values
    qs = results.groupby(
        ["data_type", "subset", "condition"]
    ).apply(lambda x: pd.Series(multipletests(x['p_value'], method="fdr_bh")[1], index=x['grna'])).reset_index().rename(columns={0: "q_value"})
    results = results.merge(qs)
    results["log_q_value"] = -np.log10(results["q_value"])
    results = results.sort_values("log_q_value", ascending=False)
    results.to_csv(os.path.join(results_dir, "{}.perturbation_assessment.csv".format(experiment)), index=False)

    g = sns.FacetGrid(data=results.sort_values("log_q_value").replace(np.inf, 360), col="data_type", row="subset", hue="condition", sharey=False, sharex=False, size=4, aspect=1)
    g.map(sns.stripplot, "gene", "log_q_value")
    g.add_legend()
    sns.despine(g.fig)
    for ax in g.axes.flat:
        ax.axhline(-np.log10(0.05))
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    g.fig.savefig(os.path.join(results_dir, "{}.perturbation_assessment.p_values.stripplot.svg".format(experiment)), bbox_inches="tight")

    # combine p-values for each gene
    cp = results.groupby(["data_type", "subset", "condition", "gene"])['p_value'].apply(lambda x: combine_pvalues(x)[1]).reset_index()
    cp["q_value"] = multipletests(cp['p_value'], method="fdr_bh")[1]
    cp['log_q_value'] = -np.log10(cp['q_value'])

    g = sns.FacetGrid(data=cp.sort_values("log_q_value").replace(np.inf, 160), col="data_type", row="subset", hue="condition", sharey=False, sharex=False, size=4, aspect=1)
    g.map(sns.stripplot, "gene", "log_q_value")
    g.add_legend()
    sns.despine(g.fig)
    for ax in g.axes.flat:
        ax.axhline(-np.log10(0.05))
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    g.fig.savefig(os.path.join(results_dir, "{}.perturbation_assessment.combined_p_values.stripplot.svg".format(experiment)), bbox_inches="tight")

    for data_type in cp["data_type"].drop_duplicates():
        cp2 = cp[(cp["data_type"] == data_type)]

        # Exemplary scatter plots
        if data_type == "CROP":
            orig_m = df
        else:
            orig_m = df_bulk

        fig, axis = plt.subplots(2, 4, figsize=(4 * 4, 2 * 4))
        axis = iter(axis.flatten())
        for subset in cp2["subset"].drop_duplicates():
            for condition in cp2["condition"].drop_duplicates():
                for f in ["tail", "head"]:
                    t = cp2[(cp2["condition"] == condition) & (cp2["subset"] == subset)].sort_values("log_q_value")
                    top_genes = getattr(t, f)(1)["gene"]
                    for gene in top_genes:
                        print(data_type, subset, condition, f, gene)
                        ax = axis.next()
                        g = orig_m[orig_m.columns[(orig_m.columns.get_level_values('condition') == condition) & (orig_m.columns.get_level_values('gene') == gene)]].mean(axis=1)
                        ctrl = orig_m[orig_m.columns[(orig_m.columns.get_level_values('condition') == condition) & (orig_m.columns.get_level_values('gene').str.contains("CTRL"))]].mean(axis=1)

                        if subset != "all_genes":
                            g = g.ix[de_genes].dropna()
                            ctrl = ctrl.ix[de_genes].dropna()

                        # x == y
                        ax.plot((0, max(g.max(), ctrl.max())), (0, max(g.max(), ctrl.max())), ls="--", lw=2, color="black", alpha=0.5)

                        # Fit lowess (for colormap)
                        l = lowess(g, ctrl, return_sorted=False)
                        dist = abs(g - l)

                        # Plot scatter
                        ax.scatter(g, ctrl, alpha=0.5, s=2, color=plt.cm.inferno(dist))
                        ax.plot((0, max(g.max(), ctrl.max())), (0, max(g.max(), ctrl.max())), ls="--", lw=2, color="black", alpha=0.5)
                        ax.set_title("; ".join([subset, condition, f, gene]))

        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "{}.perturbation_assessment.{}.top_scatters.png".format(experiment, data_type)), bbox_inches="tight", dpi=300)
        fig.savefig(os.path.join(results_dir, "{}.perturbation_assessment.{}.top_scatters.svg".format(experiment, data_type)), bbox_inches="tight")


def z_score(x):
    return (x - x.mean()) / x.std()


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


def differential_genes(df, experiment, method="pca", results_dir="results", n_genes=500):
    """
    """
    def pca(df, level="gene", filter_low=True, percentile=99):
        from sklearn.decomposition import PCA
        level_mapping = dict(zip(df.columns.names, range(len(df.columns.names))))

        if filter_low:
            df2 = df[df.sum().sort_values().tail(int(df.shape[1] * 0.25)).index]  # get top 25% covered cells
        else:
            df2 = df

        # Cells grouped by gene
        df_group = df2.T.groupby(level=[level_mapping["condition"], level_mapping[level]]).mean().T
        df_group.index = df2.index

        fitted = PCA().fit_transform(df_group)  # for genes
        r = pd.Series(fitted[:, 1], index=df2.index).sort_values()
        de_genes = r[abs(r) > np.percentile(abs(r), percentile)].index

        g = sns.clustermap(
            df2.ix[de_genes],
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
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.PCA.clustering.svg".format(experiment, n_genes)), bbox_inches="tight")

        g = sns.clustermap(
            df_group.ix[de_genes],
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
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.PCA.clustering.{}_level.svg".format(experiment, n_genes, level)), bbox_inches="tight")

        # All genes
        fig, axis = plt.subplots(1)
        sns.distplot(np.log2(1 + abs(diff)), kde=False)
        axis.set_ylabel("Genes")
        axis.set_xlabel("PC contribution")
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "{}.{}genes.PCA.all_PC_contribution.svg".format(experiment, n_genes, level)), bbox_inches="tight")

        # Significant ones
        s = diff[abs(diff) > np.percentile(abs(diff), percentile)]
        fig, axis = plt.subplots(1)
        sns.barplot(s, s.index, orient="horiz")
        axis.set_ylabel("Gene")
        axis.set_xlabel("PC contribution")
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "{}.{}genes.PCA.significant_PC_contribution.svg".format(experiment, n_genes, level)), bbox_inches="tight")

        return r

    print("Getting differential gene expression for experiment: '{}'".format(experiment))

    if method == "pca":
        diff = pca(df)
        diff.to_csv(os.path.join(results_dir, "{}.differential_expression.{}.stimutation.csv".format(experiment, method)), index=True)

    return diff


def enrich_signature(method="pca", percentile=99, results_dir="results", experiment="CROP-seq_Jurkat_TCR", n_genes=500):
    """
    """
    diff = pd.read_csv(os.path.join(results_dir, "{}.differential_expression.{}.stimutation.csv".format(experiment, method)), squeeze=True, index_col=0, header=None, names=["gene_name"])
    degs = pd.Series(diff[abs(diff) > np.percentile(abs(diff), percentile)].index)
    degs.name = "gene_name"

    enr = enrichr(degs.reset_index())
    enr.to_csv(os.path.join(results_dir, "differential_expression.{}.enrichr.csv".format(method)), index=False, encoding="utf8")

    # Plot top N terms of each library
    n = 8

    to_plot = [
        'GO_Biological_Process_2015',
        "KEGG_2016",
        "WikiPathways_2016",
        "Reactome_2016",
        "BioCarta_2016",
        "NCI-Nature_2016"]

    p = enr.ix[enr[enr['gene_set_library'].isin(to_plot)].groupby("gene_set_library")['combined_score'].nlargest(n).index.get_level_values(1)].sort_values("combined_score", ascending=False)

    fig, axis = plt.subplots(1)
    sns.barplot(data=p, y="description", x="combined_score", orient="horiz", hue="gene_set_library")
    axis.set_xlabel("Combined score")
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "differential_expression.{}.enrichr.top{}_terms.svg".format(method, n)), bbox_inches="tight")


def stimulation_signature(
        df, df_bulk, de_genes,
        experiment="CROP-seq_Jurkat_TCR", n_genes=500, method="pca", cond1="stimulated", cond2="unstimulated",
        percentile=99, results_dir="results"):
    from scipy.stats import pearsonr, spearmanr

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

    for data_type, matrix in [("CROP", df), ("Bulk", df_bulk)]:
        if data_type == "CROP":
            prefix = "{}.{}genes.{}".format(experiment, n_genes, method)
        else:
            prefix = "{}.{}.{}".format(experiment, data_type, method)

        # Cluster all cells on DE genes
        g = sns.clustermap(
            matrix.ix[de_genes],
            metric="correlation",
            z_score=0,
            vmin=-3, vmax=3,
            row_cluster=True, col_cluster=True,
            yticklabels=True, xticklabels=False,
            col_colors=get_level_colors(matrix.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.all_cells.clustering.png".format(prefix)), bbox_inches="tight", dpi=300)

        # Correlate all cells on DE genes
        g = sns.clustermap(
            matrix.ix[de_genes].corr(),
            cmap="BrBG",
            metric="correlation",
            row_cluster=True, col_cluster=True,
            yticklabels=True, xticklabels=False,
            col_colors=get_level_colors(matrix.columns),
            row_colors=get_level_colors(matrix.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.all_cells.clustering.correlation.png".format(prefix)), bbox_inches="tight", dpi=300)

        # Cluster CTRL cells on DE genes
        g = sns.clustermap(
            matrix[matrix.columns[matrix.columns.get_level_values('gene') == "CTRL"]].ix[de_genes],
            metric="correlation",
            z_score=0,
            vmin=-3, vmax=3,
            row_cluster=True, col_cluster=True,
            yticklabels=True, xticklabels=False,
            col_colors=get_level_colors(matrix.columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.control_cells.clustering.svg".format(prefix)), bbox_inches="tight")

        # Correlate CTRL cells on DE genes
        g = sns.clustermap(
            matrix[matrix.columns[matrix.columns.get_level_values('gene') == "CTRL"]].ix[de_genes].corr(),
            cmap="BrBG",
            metric="correlation",
            row_cluster=True, col_cluster=True,
            yticklabels=True, xticklabels=False,
            col_colors=get_level_colors(matrix[matrix.columns[matrix.columns.get_level_values('gene') == "CTRL"]].columns),
            row_colors=get_level_colors(matrix[matrix.columns[matrix.columns.get_level_values('gene') == "CTRL"]].columns),
            figsize=(15, 15))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.fig.savefig(os.path.join(results_dir, "{}.control_cells.clustering.correlation.svg".format(prefix)), bbox_inches="tight")

        # Cluster groups of targeted gRNAs/genes
        # get condition/gene mean expression for every gene
        for level in ["gene", "grna"]:
            matrix_level_means = matrix.T.groupby(level=["condition", level]).mean().T
            # grna level
            # cluster mean gene expression
            g = sns.clustermap(
                matrix_level_means.ix[de_genes],
                z_score=0,
                vmin=-3, vmax=3,
                row_cluster=True, col_cluster=True,
                yticklabels=False, xticklabels=True,
                col_colors=get_level_colors(matrix_level_means.columns),
                figsize=(15, 15))
            for item in g.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
                item.set_fontsize(8)
            g.fig.savefig(os.path.join(results_dir, "{}.{}_level.clustering.svg".format(prefix, level)), bbox_inches="tight")

            # order
            g = sns.clustermap(
                matrix_level_means.ix[de_genes],
                z_score=0,
                vmin=-3, vmax=3,
                row_cluster=True, col_cluster=False,
                yticklabels=False, xticklabels=True,
                col_colors=get_level_colors(matrix_level_means.columns),
                figsize=(15, 15))
            for item in g.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
                item.set_fontsize(8)
            g.fig.savefig(os.path.join(results_dir, "{}.{}_level.ordered_heatmap.svg".format(prefix, level)), bbox_inches="tight")

            # correlation
            g = sns.clustermap(
                matrix_level_means.ix[de_genes].corr(),
                cmap="BrBG",
                row_cluster=True, col_cluster=True,
                yticklabels=True, xticklabels=False,
                col_colors=get_level_colors(matrix_level_means.columns),
                row_colors=get_level_colors(matrix_level_means.columns),
                figsize=(15, 15))
            for item in g.ax_heatmap.get_yticklabels():
                item.set_rotation(0)
                item.set_fontsize(8)
            g.fig.savefig(os.path.join(results_dir, "{}.{}_level.clustering.correlation.svg".format(prefix, level)), bbox_inches="tight")

    #

    # Signature-based cell assignemnt

    # Get stimulation signature
    # 1. get mean expression of each group in signature gens
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
    for fraction in np.linspace(0, 1, 100)[1:]:
        df_sub = df * fraction
        df_sub[df_sub < 1] = 0
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

        cors.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_correlation.{}.csv".format(experiment, n_genes, fraction)))
        p_values.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_p_values.{}.csv".format(experiment, n_genes, fraction)))

    # visualize
    # sorted by max
    sigs = cors.apply(lambda x: np.argmax(x), axis=1).astype(float).sort_values()

    g = sns.clustermap(
        cors.ix[sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(cors.ix[sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_matrix.sorted.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    g = sns.clustermap(
        -np.log10(p_values).ix[sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(p_values.ix[sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.pvalue_matrix.sorted.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    for level in ['grna', 'gene']:
        c = cors.groupby(level=['condition', level]).mean()
        cs = c.apply(lambda x: np.argmax(x), axis=1).astype(float).sort_values()

        p = p_values.groupby(level=['condition', level])

        g = sns.clustermap(
            c.ix[cs.index],
            z_score=0,
            row_cluster=False, col_cluster=False,
            yticklabels=True, xticklabels=False,
            row_colors=get_level_colors(c.ix[cs.index].index))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_matrix.sorted.{}.svg".format(experiment, n_genes, level)), bbox_inches="tight")
        g.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_matrix.sorted.{}.png".format(experiment, n_genes, level)), bbox_inches="tight", dpi=300)

    # 3b. get background of signature correlations/positions
    if data_type == "CROP":
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
    i = np.random.choice(range(len(random_cors.index)), 5000)

    # sorted by max
    random_sigs = random_cors.iloc[i].apply(lambda x: np.argmax(x), axis=1).sort_values()

    g = sns.clustermap(
        random_cors.iloc[i].ix[random_sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(random_cors.iloc[i].ix[random_sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.correlation_matrix.sorted.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    g = sns.clustermap(
        -np.log10(random_p_values.iloc[i]).ix[random_sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=False, xticklabels=False,
        row_colors=get_level_colors(random_p_values.iloc[i].ix[random_sigs.index].index))
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.random_cells.pvalue_matrix.sorted.png".format(experiment, n_genes)), dpi=300, bbox_inches="tight")

    for level in ['grna', 'gene']:
        c = random_cors.groupby(level=['condition', level]).mean()
        cs = c.apply(lambda x: np.argmax(x), axis=1).sort_values()

        p = random_p_values.groupby(level=['condition', level])

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

    #

    # 3c. get signature value of each bulk sample

    # get signature matrix from bulk samples on the same genes
    df_bulk_gene_means = df_bulk.T.groupby(level=["condition", "gene"]).median().T.ix[de_genes].dropna()
    bx1 = df_bulk_gene_means[df_bulk_gene_means.columns[df_bulk_gene_means.columns.get_level_values("condition") == cond1]].median(axis=1)
    bx1.name = cond1
    bx2 = df_bulk_gene_means[df_bulk_gene_means.columns[df_bulk_gene_means.columns.get_level_values("condition") == cond2]].median(axis=1)
    bx2.name = cond2
    bulk_sign_mat = generate_signature_matrix(np.vstack([bx1, bx2]).T, n=101, bounds=(-20, 20))

    bulk_cors = list()
    bulk_p_values = list()
    for i, sample in enumerate(df_bulk.columns):
        cor, p_value = best_signature_matrix(array=df_bulk[sample].ix[df_bulk_gene_means.index], matrix=bulk_sign_mat)
        bulk_cors.append(cor)
        bulk_p_values.append(p_value)

    bulk_cors = pd.DataFrame(bulk_cors, index=df_bulk.columns)
    bulk_p_values = pd.DataFrame(bulk_p_values, index=df_bulk.columns)

    bulk_cors.to_csv(os.path.join(results_dir, "{}.{}genes.signature.Bulk.matrix_correlation.csv".format(experiment, n_genes)))
    bulk_p_values.to_csv(os.path.join(results_dir, "{}.{}genes.signature.Bulk.matrix_p_values.csv".format(experiment, n_genes)))
    bulk_cors = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.Bulk.matrix_correlation.csv".format(experiment, n_genes)), index_col=range(4))
    bulk_p_values = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.Bulk.matrix_p_values.csv".format(experiment, n_genes)), index_col=range(4))

    # visualize
    # sorted by max
    sigs = bulk_cors.apply(lambda x: np.argmax(x), axis=1).astype(float).sort_values()

    g = sns.clustermap(
        bulk_cors.ix[sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=True, xticklabels=False,
        row_colors=get_level_colors(bulk_cors.ix[sigs.index].index))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
        item.set_fontsize(8)
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.Bulk.correlation_matrix.sorted.svg".format(experiment, n_genes)), bbox_inches="tight")

    g = sns.clustermap(
        -np.log10(bulk_p_values).ix[sigs.index],
        z_score=0,
        row_cluster=False, col_cluster=False,
        yticklabels=True, xticklabels=False,
        row_colors=get_level_colors(bulk_p_values.ix[sigs.index].index))
    for item in g.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
        item.set_fontsize(8)
    g.savefig(os.path.join(results_dir, "{}.{}genes.signature.Bulk.pvalue_matrix.sorted.svg".format(experiment, n_genes)), bbox_inches="tight")

    for level in ['grna', 'gene']:
        c = bulk_cors.groupby(level=['condition', level]).mean()
        cs = c.apply(lambda x: np.argmax(x), axis=1).astype(float).sort_values()

        p = bulk_p_values.groupby(level=['condition', level])

        g = sns.clustermap(
            c.ix[cs.index],
            z_score=0,
            row_cluster=False, col_cluster=False,
            yticklabels=True, xticklabels=False,
            row_colors=get_level_colors(c.ix[cs.index].index))
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
            item.set_fontsize(8)
        g.savefig(os.path.join(results_dir, "{}.{}genes.signature.Bulk.correlation_matrix.sorted.{}.svg".format(experiment, n_genes, level)), bbox_inches="tight")

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
    for data_type, matrix in [("CROP", df), ("Bulk", df_bulk)]:
        if data_type == "CROP":
            index = ['condition', 'replicate', 'cell', 'grna', 'gene']
            cors = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.matrix_correlation.csv".format(experiment, n_genes, "all_cells")))
            cors['replicate'] = cors['replicate'].astype(str)
            cors = cors.set_index(index)
            p_values = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.matrix_p_values.csv".format(experiment, n_genes, "all_cells")))
            p_values['replicate'] = p_values['replicate'].astype(str)
            p_values = p_values.set_index(index)
        else:
            index = ['condition', 'sample_name', 'grna', 'gene']
            cors = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.matrix_correlation.csv".format(experiment, n_genes, data_type))).set_index(index)
            p_values = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.matrix_p_values.csv".format(experiment, n_genes, data_type))).set_index(index)

        sigs = cors.apply(lambda x: np.argmax(x), axis=1)
        sigs.name = "signature"

        # get "uncorrelated" fraction (noise)
        res = 1 - cors.max(axis=1)
        res.name = "residual"

        # reads per cell
        sigs = pd.merge(sigs.reset_index(), res.reset_index()).set_index(index)

        # reads per cell
        if data_type == "CROP":
            s = df.sum(axis=0)
        elif data_type == "Bulk":
            s = df_bulk.sum(axis=0)
        s.name = "reads_per_cell"
        sigs = pd.merge(sigs.reset_index(), s.reset_index()).set_index(index)

        # get minimum correlation
        mean = np.mean
        for metric in ["min", "mean", "max"]:
            m = cors.apply(eval(metric), axis=1)
            m.name = "{}_corr".format(metric)
            sigs = pd.merge(sigs.reset_index(), m.reset_index()).set_index(index)

        # save
        sigs.to_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.all_cells.correlation.csv".format(experiment, n_genes, data_type)))
        sigs = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.all_cells.correlation.csv".format(experiment, n_genes, data_type))).set_index(index)

        # Plot heatmaps of metrics per cell sorted by signature (to match the correlation matrices per cell)
        for metric in ["residual", "min_corr", "mean_corr", "max_corr"]:
            print(metric)
            fig, axis = plt.subplots(1, figsize=(2, 12))
            sns.heatmap(sigs.sort_values("signature")[[metric]], ax=axis, cmap="inferno", xticklabels=False, yticklabels=False, vmin=0, vmax=1)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.all_cells.{}.heatmap_strip.png".format(experiment, n_genes, data_type, metric)), bbox_inches="tight", dpi=300)
        metric = "reads_per_cell"
        fig, axis = plt.subplots(1, figsize=(2, 12))
        sns.heatmap(np.log10(sigs.sort_values("signature")[[metric]]), ax=axis, cmap="inferno", xticklabels=False, yticklabels=False)
        fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.all_cells.{}.heatmap_strip.png".format(experiment, n_genes, data_type, metric)), bbox_inches="tight", dpi=300)

        # pairwise variable distribution
        g = sns.pairplot(sigs.reset_index(), vars=sigs.columns, hue="condition")
        g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.{}.all_metrics.pairwise_distributions.png".format(experiment, n_genes, data_type)), bbox_inches="tight", dpi=300)

        # annotate KO genes with signature at both levels
        for level in ["gene", "grna"]:
            sigs_mean = sigs.astype(float).groupby(level=['condition', level]).mean()
            sigs_mean["n_cells"] = sigs.groupby(level=['condition', level]).apply(len)

            sigs_mean['signature'] = sigs_mean['signature'].replace(0, 1)  # to avoid -inf fold changes, could be done differently too

            # Plot heatmaps of metrics per cell sorted by signature (to match the correlation matrices per cell)
            for metric in ["residual", "min_corr", "mean_corr", "max_corr"]:
                print(metric)
                fig, axis = plt.subplots(1, figsize=(2, 12))
                sns.heatmap(sigs_mean.sort_values("signature")[[metric]], ax=axis, cmap="inferno", xticklabels=False, yticklabels=False, vmin=0, vmax=1)
                fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.{}.{}.heatmap_strip.png".format(experiment, n_genes, data_type, level, metric)), bbox_inches="tight", dpi=300)
            metric = "reads_per_cell"
            fig, axis = plt.subplots(1, figsize=(2, 12))
            sns.heatmap(np.log10(sigs_mean.sort_values("signature")[[metric]]), ax=axis, cmap="inferno", xticklabels=False, yticklabels=False)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.{}.{}.heatmap_strip.png".format(experiment, n_genes, data_type, level, metric)), bbox_inches="tight", dpi=300)

            # add distance from CTRL
            for cond in [cond1, cond2]:
                ko = sigs_mean.loc[sigs_mean.index.get_level_values('condition') == cond, "signature"]
                ctrl = sigs_mean.loc[(sigs_mean.index.get_level_values('condition') == cond) & (sigs_mean.index.get_level_values(level).str.contains("CTRL")), "signature"].mean()
                sigs_mean.loc[sigs_mean.index.get_level_values('condition') == cond, "abs_change"] = ko - ctrl
                sigs_mean.loc[sigs_mean.index.get_level_values('condition') == cond, "log_fold_change"] = np.log2(ko / ctrl)

            # save
            sigs_mean = sigs_mean.sort_index(level='condition')
            sigs_mean = sigs_mean.sort_values(['signature'])
            sigs_mean.to_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.{}_means.annotated.csv".format(experiment, n_genes, data_type, level)), index=True)
            sigs_mean = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.{}_means.annotated.csv".format(experiment, n_genes, data_type, level)), index_col=[0, 1])

            #

            # Filter out genes with less than n of cells
            if data_type == "CROP":
                sigs_mean = sigs_mean[sigs_mean["n_cells"] >= 10]

            fig, axis = plt.subplots(1, figsize=(10, 8))
            sns.stripplot(x=sigs_mean['signature'], y=sigs_mean.index, orient="horiz", ax=axis)
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.all_cells.mean_{}_signature.strip.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            # Given the known condition of each group, what is the deviation from that?
            # plot as rank of mean
            fig, axis = plt.subplots(1, figsize=(10, 8))
            axis.scatter(sigs_mean['signature'].rank(ascending=False), sigs_mean['signature'])
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.all_cells.mean_{}_signature.rank.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

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
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.all_cells.mean_{}_signature.violinplot.sorted_st.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            fig, axis = plt.subplots(1, figsize=(16, 8))
            sns.violinplot(
                x="condition", y="signature", hue="gene",
                data=p,
                cut=0,
                order=[cond1, cond2],
                hue_order=sigs_mean.reset_index().sort_values(["signature"])["gene"].drop_duplicates(),
                ax=axis)
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.all_cells.mean_{}_signature.violinplot.sorted_un.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            fig, axis = plt.subplots(1, figsize=(16, 8))
            sns.violinplot(x="gene", y="signature", hue="condition", cut=0, data=sigs.reset_index(), ax=axis)
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.all_cells.mean_{}_signature.violinplot2.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            # Make heatmap sorted by median signature position per knockout

            # Group by stimulation / gRNA, get mean expression
            level_matrix = matrix.ix[de_genes].dropna().T.groupby(level=['condition', level]).mean()

            # cluster
            g = sns.clustermap(
                level_matrix.T, z_score=0,
                col_colors=get_level_colors(level_matrix.index),
                metric='correlation',
                row_cluster=True, col_cluster=True,
                xticklabels=True, yticklabels=True,
                figsize=(15, 15))
            for item in g.ax_heatmap.get_yticklabels():
                item.set_rotation(0)
            for item in g.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
            g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.{}.{}_group_means.clustered.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            # cluster groups, sort genes
            g = sns.clustermap(
                level_matrix.ix[sigs_mean.sort_values("signature").index].T, z_score=0,
                col_colors=get_level_colors(level_matrix.ix[sigs_mean.sort_values("signature").index].index),
                metric='correlation',
                row_cluster=False, col_cluster=True,
                xticklabels=True, yticklabels=True,
                figsize=(15, 15))
            for item in g.ax_heatmap.get_yticklabels():
                item.set_rotation(0)
            for item in g.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
            g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.{}.{}_group_means.sorted_genes_clustered_groups.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            # sort by signature
            clust = sns.clustermap(
                level_matrix.ix[sigs_mean.sort_values("signature").index], z_score=1,
                row_colors=get_level_colors(sigs_mean.sort_values("signature").index),
                metric='correlation',
                row_cluster=False, col_cluster=True,
                xticklabels=True, yticklabels=True,
                figsize=(6, 8.62948158106742))
            for item in clust.ax_heatmap.get_yticklabels():
                item.set_rotation(0)
            for item in clust.ax_heatmap.get_xticklabels():
                item.set_rotation(90)
            clust.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.{}.{}_group_means.sorted_signature.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            # Strip plot of number of cells
            p = sigs_mean.sort_values("signature")

            fig, axis = plt.subplots(1, figsize=(8, 8))
            axis.scatter([1] * p.shape[0], range(p.shape[0]), s=p['n_cells'])
            [axis.text(1.0005, k, s=str(int(x))) for k, x in enumerate(p['n_cells'].values)]
            [axis.text(1.01, k, s=str(x)) for k, x in enumerate(p.index.get_level_values('condition'))]
            [axis.text(1.02, k, s=str(x)) for k, x in enumerate(p.index.get_level_values(level))]
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.{}_level.cells_per_group.bubbles.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            # Stripe with relative change
            fig, axis = plt.subplots(1, figsize=(8, 8))
            sns.heatmap(p[["log_fold_change"]], ax=axis)
            fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.{}_level.cells_per_group.stripe.svg".format(experiment, n_genes, data_type, level)), bbox_inches="tight")

            # Barplots of difference compared with CTRL
            for metric in ["abs_change", "log_fold_change"]:
                fig, axis = plt.subplots(1, figsize=(4, 4))
                sns.barplot(sigs_mean[metric], sigs_mean[metric].index, orient="horiz", order=sigs_mean[metric].index, ax=axis)
                axis.set_xlabel("Signature deviation from control (propensity to cause stimulation)")
                sns.despine(fig)
                fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.{}.{}_mean_group_signature_deviation.{}.barplot.svg".format(experiment, n_genes, data_type, level, metric)), bbox_inches="tight")

            if data_type == "CROP":
                # Single-cell matrix sorted in same way as above
                signature_dict = sigs_mean.sort_values("signature")['signature'].to_dict()

                df2 = df.copy().T
                df2.index = df2.index.droplevel(['cell', 'replicate', 'grna'])
                df2['sortby'] = df2.index.to_series().map(signature_dict).tolist()
                df2 = df2.reset_index().sort_values(["sortby", "condition", level]).drop(['sortby'], axis=1).set_index(['condition', level]).T

                # sort by signature
                g = sns.clustermap(
                    df2.ix[de_genes].T,
                    z_score=1, vmin=-1.5, vmax=1.5,
                    row_colors=get_level_colors(df2.columns),
                    metric='correlation',
                    row_cluster=False, col_cluster=True,
                    xticklabels=True, yticklabels=True,
                    figsize=(15, 15))
                for item in g.ax_heatmap.get_yticklabels():
                    item.set_rotation(0)
                for item in g.ax_heatmap.get_xticklabels():
                    item.set_rotation(90)
                g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.{}.single_cells.sorted_signature.png".format(experiment, n_genes, data_type)), bbox_inches="tight", dpi=300)
                # g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.{}.single_cells.sorted_signature.svg".format(experiment, n_genes, data_type)), bbox_inches="tight")

    #

    # Compare CROP and Bulk RNA-seq head-to-head at the two levels
    for level in ["gene", "grna"]:
        sigs_mean = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.{}_means.annotated.csv".format(experiment, n_genes, "CROP", level)))
        sigs_mean['data_type'] = "CROP"
        bulk_sigs_mean = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.{}_means.annotated.csv".format(experiment, n_genes, "Bulk", level)))
        bulk_sigs_mean['data_type'] = "Bulk"

        t = sigs_mean.append(bulk_sigs_mean)

        if level == 'grna':
            t['grna'] = t['grna'].str.replace("Tcrlibrary_", "")

        # Scatter plot of signature positioning
        # Scatter with fold changes relative to CTRL
        condition_dict = dict(zip(t['condition'].drop_duplicates(), sns.color_palette("colorblind")))

        fig0, axis0 = plt.subplots(1, 2, figsize=(5 * 2, 5 * 1))
        fig1, axis1 = plt.subplots(1, 2, figsize=(5 * 2, 5 * 1))
        t_pivot = pd.pivot_table(t, index=['condition', level], columns="data_type", values="signature").dropna()
        p = pearsonr(t_pivot["CROP"], t_pivot["Bulk"])[0]
        s = spearmanr(t_pivot["CROP"], t_pivot["Bulk"])[0]

        axis0[0].scatter(t_pivot["CROP"], t_pivot["Bulk"], color=[condition_dict[k] for k in t_pivot.index.get_level_values("condition")], alpha=0.8)
        axis0[0].text(t_pivot["CROP"].max(), t_pivot["Bulk"].min(), "Pearson: {0:0.3f}\nSpearman: {1:0.3f}".format(p, s))
        axis1[0].scatter(t_pivot["CROP"], t_pivot["Bulk"], color=[condition_dict[k] for k in t_pivot.index.get_level_values("condition")], alpha=0.8)
        axis1[0].text(t_pivot["CROP"].max(), t_pivot["Bulk"].min(), "Pearson: {0:0.3f}\nSpearman: {1:0.3f}".format(p, s))
        for g in t_pivot.index:
            axis1[0].text(t_pivot.ix[g]["CROP"], t_pivot.ix[g]["Bulk"], g[1], color=condition_dict[g[0]], alpha=0.8)

        t_pivot = pd.pivot_table(t, index=['condition', level], columns="data_type", values="log_fold_change").dropna()
        p = pearsonr(t_pivot["CROP"], t_pivot["Bulk"])[0]
        s = spearmanr(t_pivot["CROP"], t_pivot["Bulk"])[0]

        axis0[1].scatter(t_pivot["CROP"], t_pivot["Bulk"], color=[condition_dict[k] for k in t_pivot.index.get_level_values("condition")], alpha=0.8)
        axis0[1].text(t_pivot["CROP"].max(), t_pivot["Bulk"].min(), "Pearson: {0:0.3f}\nSpearman: {1:0.3f}".format(p, s))
        axis1[1].scatter(t_pivot["CROP"], t_pivot["Bulk"], color=[condition_dict[k] for k in t_pivot.index.get_level_values("condition")], alpha=0.8)
        axis1[1].text(t_pivot["CROP"].max(), t_pivot["Bulk"].min(), "Pearson: {0:0.3f}\nSpearman: {1:0.3f}".format(p, s))
        for g in t_pivot.index:
            axis1[1].text(t_pivot.ix[g]["CROP"], t_pivot.ix[g]["Bulk"], g[1], color=condition_dict[g[0]], alpha=0.8)

        axis0[0].set_title("Predicted signature position")
        axis0[1].set_title("Log fold change of TCR pathway influence")
        axis1[0].set_title("Predicted signature position")
        axis1[1].set_title("Log fold change of TCR pathway influence")
        for q in [axis0, axis1]:
            for ax in q:
                ax.set_xlabel("CROP-seq")
                ax.set_ylabel("RNA-seq on bulk population of KOs")
        sns.despine(fig0)
        fig0.savefig(os.path.join(results_dir, "{}.predicted_signature.CROP_vs_Bulk.{}.scatter.svg".format(experiment, level)), bbox_inches="tight")
        sns.despine(fig1)
        fig1.savefig(os.path.join(results_dir, "{}.predicted_signature.CROP_vs_Bulk.{}.scatter+text.svg".format(experiment, level)), bbox_inches="tight")

    #

    # Rarefaction analysis:
    # we're working with the CROP-seq data for each single-cell:

    # The assumption is that the quality of the transcriptome of each single-cell would be the same as currently
    # and therefore we'll vary only the number of discovered cells per gRNA or gene.
    # For each subsample we'll compare the CROP-seq data with the bulk and observe several metrics
    index = ['condition', 'sample_name', 'grna', 'gene']
    bulk_sigs = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.{}.all_cells.correlation.csv".format(experiment, n_genes, "Bulk"))).set_index(index)

    # Start subsampling in 100 fractions of the data
    n_iter = 100
    rare_metrics = pd.DataFrame()
    for i in range(n_iter):
        for umi_fraction in np.linspace(0.0, 1.0, 100)[1:]:
            data_type = "CROP"
            index = ['condition', 'replicate', 'cell', 'grna', 'gene']
            cors = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_correlation.{}.csv".format(experiment, n_genes, umi_fraction))).set_index(index)
            sigs = cors.apply(lambda x: np.argmax(x), axis=1)
            sigs.name = "signature"
            sigs = pd.DataFrame(sigs)
            for condition in [cond1]:
                # Compare aggregate signatures at both levels:
                for level in ["gene", "grna"]:

                    bulk_sigs_mean = bulk_sigs.astype(float).groupby(level=['condition', level]).mean()
                    bulk_sigs_mean['signature'] = bulk_sigs_mean['signature'].replace(0, 1)  # to avoid -inf fold changes, could be done differently too
                    bulk_sigs_mean['data_type'] = "Bulk"

                    for cell_fraction in np.linspace(0.0, 1.0, 100)[1:]:
                        print(i, umi_fraction, cell_fraction, condition, level)
                        total = sigs.shape[0]
                        # Subsample a fraction of all cells (maintains the distribution of recovered cells)
                        chosen_indexes = np.random.choice(sigs.index, int(np.round(cell_fraction * total)))

                        sigs_mean = sigs.ix[chosen_indexes].astype(float).groupby(level=['condition', level]).mean()
                        sigs_mean['signature'] = sigs_mean['signature'].replace(0, 1)  # to avoid -inf fold changes, could be done differently too
                        sigs_mean['data_type'] = "CROP"

                        # Match CROP and Bulk
                        t = sigs_mean.append(bulk_sigs_mean).reset_index()
                        if level == 'grna':
                            t['grna'] = t['grna'].str.replace("Tcrlibrary_", "")

                        # Measure correlations
                        t_pivot = pd.pivot_table(t, index=['condition', level], columns="data_type", values="signature").dropna()
                        p = pearsonr(t_pivot["CROP"], t_pivot["Bulk"])[0]

                        rare_metrics = rare_metrics.append(pd.Series(
                            [i, condition, level, umi_fraction, cell_fraction, p],
                            index=["iteration", "condition", "level", "umi_fraction", "cell_fraction", "pearson"]), ignore_index=True)

        rare_metrics.to_csv(os.path.join(results_dir, "{}.signature.rarefaction_analysis.{}_iterations.csv".format(experiment, i)), index=False)

    rare_metrics = pd.read_csv(os.path.join(results_dir, "{}.signature.rarefaction_analysis.{}_iterations.csv".format(experiment, n_iter)))

    # Plot metrics

    # Only cell number subsampling
    rare_metrics2 = rare_metrics[rare_metrics["umi_fraction"] == rare_metrics["umi_fraction"].max()].drop(["iteration", "umi_fraction"], axis=1)

    # cell states combined CI
    g = sns.factorplot(
        x="cell_fraction", y="pearson", row="level", hue="condition", data=rare_metrics2,
        ci=.95, palette="colorblind", s=1, alpha=0.8, figsize=(20, 10))
    for ax in g.axes.flat:
        ax.set_xlabel("Fraction of cells sampled")
        ax.set_ylabel("Correlation with RNA-seq of Bulk KO cell line")
        ax.set_ylim(0.5, 1.0),
        ax.set_xticklabels(ax.get_xticklabels(), rotation="90")
    g.savefig(os.path.join(results_dir, "{}.signature.rarefaction_analysis.{}_iterations.cell_sampling.combined_ci.svg".format(experiment, n_iter)), bbox_inches="tight")

    # Only UMI subsampling
    rare_metrics2 = rare_metrics[rare_metrics["cell_fraction"] == rare_metrics["cell_fraction"].max()].drop(["iteration", "cell_fraction"], axis=1)

    # cell states combined CI
    g = sns.factorplot(
        x="umi_fraction", y="pearson", row="level", hue="condition", data=rare_metrics2,
        ci=.95, palette="colorblind", s=1, alpha=0.8, figsize=(20, 10))
    for ax in g.axes.flat:
        ax.set_xlabel("Fraction of cells sampled")
        ax.set_ylabel("Correlation with RNA-seq of Bulk KO cell line")
        # ax.set_ylim(0.5, 1.0),
        ax.set_xticklabels(ax.get_xticklabels(), rotation="90")
    g.savefig(os.path.join(results_dir, "{}.signature.rarefaction_analysis.{}_iterations.umi_sampling.combined_ci.svg".format(experiment, n_iter)), bbox_inches="tight")

    # Cell number and read subsampling
    cmaps = ["viridis", "inferno", "Spectral_r", "GnBu", "Blues", "PuBuGn", "YlGnBu", "Greys"]
    cmaps = ["Spectral_r", "YlGnBu"]
    for cmap in cmaps:
        fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 3 * 1))
        for i, level in enumerate(rare_metrics["level"].drop_duplicates()):
            print(cmap, i)
            # make pivot
            rare_pivot = pd.pivot_table(
                rare_metrics[(rare_metrics['level'] == level) & (rare_metrics['condition'] == condition)],
                index="cell_fraction", columns="umi_fraction", values="pearson").sort_index(ascending=False)

            from scipy import ndimage
            data = ndimage.gaussian_filter(rare_pivot.fillna(rare_pivot.min().min()), sigma=3.0, order=0)

            im1 = axis[i].imshow(
                data,
                cmap=cmap, vmin=0.0, vmax=1.0)
            im2 = axis[i].contour(data, levels=np.linspace(0.5, 1.0, 11))  # np.logspace(-1, 0, 10, base=2)
            for line in im2.ax.get_children():
                try:
                    line.set_color("black")
                    line.set_linewidth(0.5)
                except:
                    continue
            # axis[i].clabel(im2, inline=1, fontsize=10, colors="black", fmt="%.2f")
            axis[i].set_xlim((8, 98))
            axis[i].set_ylim((98, 0))

        cbar = plt.colorbar(im1)
        cbar.ax.set_ylabel('correlation to \n bulk RNA-seq')

        sns.despine(fig, top=True, right=True, left=True, bottom=True)
        fig.savefig(os.path.join(
            results_dir, "..",
            "{}.signature.rarefaction_analysis.{}_iterations.{}.joint_sampling.heatmap.{}.png".format(
                experiment, n_iter, condition, cmap)), bbox_inches="tight", dpi=300)
        fig.savefig(os.path.join(
            results_dir, "..",
            "{}.signature.rarefaction_analysis.{}_iterations.{}.joint_sampling.heatmap.{}.svg".format(
                experiment, n_iter, condition, cmap)), bbox_inches="tight")


def intra_variability(df, df_bulk, de_genes, experiment="CROP-seq_Jurkat_TCR", results_dir="results"):
    """
    Measure between-gRNA, intra-gene variability
    """
    from sklearn.metrics.pairwise import pairwise_distances
    from scipy.stats import mannwhitneyu
    from scipy.stats import norm
    from scipy.stats import combine_pvalues
    from statsmodels.sandbox.stats.multicomp import multipletests

    def vertical_mean_line(x, **kwargs):
        plt.axvline(x.mean(), **kwargs)

    def z_score(x):
        return (x - x.mean()) / x.std()

    # Pairwise distances between gRNAs
    results = pd.DataFrame()
    for data_type, tmp_df in [("CROP", df), ("Bulk", df_bulk)]:
        for subset, index in [("all_genes", df.index), ("sig_genes", de_genes)]:

            # groupby gRNA, reduce to median
            df_guide = tmp_df.ix[index].dropna().T.groupby(level=["condition", "gene", "grna"]).mean().T

            # compute all pairwise distances
            # euc = pd.DataFrame(np.triu(squareform(pdist(df_guide.T, metric='euclidean'))), index=df_guide.columns, columns=df_guide.columns).replace(0, np.nan)
            euc = pd.DataFrame(np.triu(pairwise_distances(df_guide.T, metric='l2')), index=df_guide.columns, columns=df_guide.columns).replace(0, np.nan)

            # ignore self-distances
            np.fill_diagonal(euc.values, pd.np.nan)

            # get intra-gene and between-genes gRNA distances
            for condition in euc.columns.levels[0]:
                for gene in euc.columns.levels[1]:
                    for mask_name, mask_y in [  # either intra-gene or inter-gene distances
                            ("intra", (euc.index.get_level_values('gene') == gene) & (euc.index.get_level_values('condition') == condition)),
                            ("inter", (euc.index.get_level_values('gene') != gene) & (euc.index.get_level_values('condition') == condition))]:
                        # get values
                        series = pd.DataFrame(euc.loc[
                            (euc.index.get_level_values('gene') == gene) & (euc.index.get_level_values('condition') == condition),
                            mask_y
                        ].values.flatten(), columns=["distances"]).dropna()
                        # label data
                        series["data_type"] = data_type
                        series["subset"] = subset
                        series["condition"] = condition
                        series["gene"] = gene
                        series["relation"] = mask_name
                        results = results.append(series)

    results.to_csv(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.euc.csv".format(experiment)), index=False)
    results = pd.read_csv(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.euc.csv".format(experiment)))

    # Test differences between type of gRNA pairs
    fig, axis = plt.subplots(1)
    results.groupby(["data_type", "subset", "condition"]).apply(
        lambda x: -np.log10(mannwhitneyu(x[x["relation"] == "intra"]["distances"], x[x["relation"] == "inter"]["distances"])[1])
    ).sort_values().plot(kind='bar', ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.euc.p_values.svg".format(experiment)), bbox_inches="tight")

    # Test differences between data types
    fig, axis = plt.subplots(1)
    results.groupby(["relation", "subset", "condition"]).apply(
        lambda x: -np.log10(mannwhitneyu(x[x["data_type"] == "CROP"]["distances"], x[x["data_type"] == "Bulk"]["distances"])[1])
    ).sort_values().plot(kind='bar', ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.euc.data_type.p_values.svg".format(experiment)), bbox_inches="tight")

    # Plot distributions
    for sub_name in results["subset"].drop_duplicates():
        results2 = results[results["subset"] == sub_name]

        # For all together
        g = sns.FacetGrid(data=results2, row="data_type", col="condition", hue="relation")
        g.map(sns.distplot, "distances", kde=True, bins=200, hist=False, kde_kws={"shade": True})
        g.map(vertical_mean_line, "distances")
        g.add_legend()
        sns.despine(g.fig)
        g.fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.euc.distplot.svg".format(experiment, sub_name)), bbox_inches="tight")

        # barplot
        g = sns.FacetGrid(data=results2, row="data_type", col="condition")
        g.map(sns.barplot, "relation", "distances")
        g.map(vertical_mean_line, "distances")
        g.add_legend()
        sns.despine(g.fig)
        g.fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.euc.barplot.svg".format(experiment, sub_name)), bbox_inches="tight")

        # barplot (CROP vs Bulk)
        fig, axis = plt.subplots(1)
        sns.barplot(data=results2, x="condition", y="distances", hue="data_type")
        sns.despine(fig)
        fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.euc.barplot.svg".format(experiment, sub_name)), bbox_inches="tight")

        g = sns.FacetGrid(data=results2, row="data_type", col="condition")
        g.map(sns.violinplot, "relation", "distances")
        sns.despine(g.fig)
        g.fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.euc.violinplot.svg".format(experiment, sub_name)), bbox_inches="tight")

        # For each data type separately
        for data_type in results2['data_type'].drop_duplicates():
            results3 = results2.loc[(results2["relation"] == "intra") & (results2["data_type"] == data_type), :]

            # For each gene
            g = sns.FacetGrid(data=results3, col="gene", col_wrap=5, hue="condition", sharey=False, xlim=(0, 70))
            g.map(sns.distplot, "distances", kde=True, hist=False, kde_kws={"shade": True})
            g.map(vertical_mean_line, "distances")
            g.map(sns.rugplot, "distances")
            g.add_legend()
            sns.despine(g.fig)
            g.fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.{}.euc.distplot.gene_level.svg".format(experiment, sub_name, data_type)), bbox_inches="tight")

            # Call significantly variable genes (between gRNAs)
            # fit gaussian on pairwise distances between control gRNAs (null)
            params = norm.fit(z_score(results3.loc[(results3["gene"] == "CTRL"), 'distances']))
            # z-score pairwise distances within all genes for each gene
            z_scores = z_score(results3.loc[(results3["relation"] == "intra", 'distances')])
            # get p-values
            results3.loc[:, 'p_value'] = norm.sf(abs(z_scores), *params) * 2  # twosided p-value
            results3.loc[:, 'q_value'] = multipletests(results3.loc[:, 'p_value'], method="fdr_bh")[1]
            results3.to_csv(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.{}.euc.metrics.csv".format(experiment, sub_name, data_type)), index=False)

            # combine p-values for each gene
            diff = results3.groupby(["condition", "gene"])['q_value'].apply(lambda x: combine_pvalues(x)[1]).reset_index()
            diff['log_q_value'] = -np.log10(diff['q_value'])

            # get variability per gene
            dist_mean = results3.groupby(["condition", "gene"])['distances'].mean()
            dist_mean.name = "distance_mean"
            dist_std = results3.groupby(["condition", "gene"])['distances'].std()
            dist_std.name = "distance_std"
            dist = pd.DataFrame([dist_mean, dist_std]).T

            # get gRNA efficiency and cell number
            annot = pd.read_csv(os.path.join("metadata", "guide_annotation.csv"))
            grna_var = annot.groupby(["gene"])["specificity_score", "efficiency_score"].std().reset_index()

            # add cell number
            n_cells = df.T.groupby(level=["condition", "gene", "grna"]).apply(len).groupby(level=['condition', 'gene']).std()
            n_cells.name = "n_cells"

            # put all together
            diff = pd.merge(
                pd.merge(
                    pd.merge(diff, dist.reset_index(), on=["condition", "gene"]),
                    grna_var, on="gene"),
                n_cells.reset_index(), on=["condition", "gene"]).set_index(["condition", "gene"])

            # Rank vs distance with p-value as size
            cmap = dict(zip(set(diff.index.get_level_values('condition')), sns.color_palette("colorblind")))

            fig, axis = plt.subplots(1)
            axis.scatter(diff["distance_mean"].rank(), diff["distance_mean"], s=20 + (5 * diff["log_q_value"]), color=[cmap[x] for x in diff.index.get_level_values('condition')])
            for i in diff.index:
                axis.text(diff["distance_mean"].rank().ix[i], diff.ix[i]["distance_mean"], " ".join(i))
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.{}.euc.metrics.rank.svg".format(experiment, sub_name, data_type)), bbox_inches="tight")

            fig, axis = plt.subplots(1)
            sns.stripplot(x="gene", y="distances", hue="condition", data=results3.sort_values("distances"), ax=axis)
            axis.axhline(np.percentile(results3[results3['gene'] == "CTRL"]['distances'], 99))
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.{}.euc.metrics.stripplot.svg".format(experiment, sub_name, data_type)), bbox_inches="tight")

            # "Volcano plot"
            fig, axis = plt.subplots(1)
            axis.scatter(diff["distance_mean"], diff["log_q_value"])
            for i in diff.index:
                axis.text(diff.ix[i]["distance_mean"], diff.ix[i]["log_q_value"], " ".join(i))
            sns.despine(fig)
            fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.{}.euc.metrics.volcano.svg".format(experiment, sub_name, data_type)), bbox_inches="tight")

            # Vizualize all pairwise relationships
            if data_type == "Bulk":
                diff = diff.drop("distance_std", axis=1)
            g = sns.pairplot(diff.dropna().reset_index(), hue="condition")
            g.fig.savefig(os.path.join(results_dir, "{}.intra_gene_variation_grna_level.{}.{}.euc.metrics.pairplot.svg".format(experiment, sub_name, data_type)), bbox_inches="tight")


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


def inspect_bulk(df, df_bulk, de_genes, de_genes_bulk):
    """
    """
    quant_types = [("bitseq", df_bulk)]

    for quant_type, exp_matrix in quant_types:
        print(quant_type)

        # Boxplots of expression
        fig, axis = plt.subplots(1)
        sns.boxplot(data=pd.melt(exp_matrix), x="grna", y="value", hue="condition", ax=axis)
        fig.savefig(os.path.join("results", "bulk", "bulk_samples.qc.{}.expression_boxplots.png".format(quant_type)), dpi=300, bbox_inches="tight")

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


def compare_bulk(df, df_bulk, de_genes, experiment="CROP-seq_Jurkat_TCR", cond1="stimulated", cond2="unstimulated"):
    """
    Compare Bulk RNA-seq data with single-cell
    """
    from scipy.stats import pearsonr, spearmanr
    from sklearn.metrics import mean_squared_error
    from math import sqrt

    def rmse(x, y):
        return sqrt(mean_squared_error(x, y)), pd.np.nan

    def n_cells(x):
        return x.shape[1], pd.np.nan

    # Tasks:

    # -1) Correlate Bulk and grna/gene level aggregated
    # 0). Compare expression in CTRL cells/samples only
    # 1). Compare expression levels across all genes
    # 2). Compare expression levels in signature genes
    # 2b). Compare expression in selected markers
    # 3). Compare signature (recall, agreement)
    # 4). Compare differential genes in KO vs CTRL (recall, agreement, enrichemnts)
    # x). Each single-cell vs Bulk
    single_quant_types = [("sc", df)]

    bulk_quant_types = [("bitseq", df_bulk)]

    # 0).
    for i, condition in enumerate(df.columns.get_level_values('condition').drop_duplicates()):
        fig, axis = plt.subplots(5, 5, figsize=(5 * 3, 5 * 3))
        axis = iter(axis.flatten())
        for gene_group, gene_filter in [("all_genes", df.index), ("de_genes", de_genes)]:
            for single_quant_type, df in single_quant_types:
                for bulk_quant_type, df_bulk in bulk_quant_types:
                    # remove gRNAs
                    df = df[~df.index.str.contains("library|CTRL")]

                    # remove ribosomal, mitochondrial genes
                    df = df[(~df.index.str.contains("^RP.*")) & (~df.index.str.contains("^MT-"))]
                    df_bulk = df_bulk[(~df_bulk.index.str.contains("^RP.*")) & (~df_bulk.index.str.contains("^MT-"))]

                    # align indices
                    df = df.ix[df_bulk.index].dropna()
                    df_bulk = df_bulk.ix[df.index].dropna()
                    de_genes = [x for x in de_genes if x in df_bulk.index]

                    a = df[df.columns[
                        (df.columns.get_level_values('condition') == condition) &
                        (df.columns.get_level_values('gene') == "CTRL")]].median(axis=1).ix[gene_filter].dropna()

                    b = df_bulk[df_bulk.columns[
                        (df_bulk.columns.get_level_values('condition') == condition) &
                        (df_bulk.columns.get_level_values('gene') == "CTRL")]].median(axis=1).ix[gene_filter].dropna()

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
    # Run
    comparisons = pd.DataFrame()
    for gene_group, gene_filter in [("all_genes", df_bulk.index)]:
        # At gene level or at grna level
        for level in ["gene", "grna"]:
            # fig, axis = plt.subplots(
            #     len(df.columns.get_level_values('condition').drop_duplicates()),
            #     len(df.columns.get_level_values(level).drop_duplicates()),
            #     figsize=(
            #         4 * len(df.columns.get_level_values(level).drop_duplicates()),
            #         4 * len(df.columns.get_level_values('condition').drop_duplicates())
            #     ), sharex=True, sharey=True
            # )
            q = pd.DataFrame()
            for i, condition in enumerate(df.columns.get_level_values('condition').drop_duplicates()):
                # Compare within each knockout
                for j, ko in enumerate(df.columns.get_level_values(level).drop_duplicates()):
                    print(gene_group, level, condition, ko)

                    s = df.columns[(df.columns.get_level_values("condition") == condition) & (df.columns.get_level_values(level) == ko)]
                    b = df_bulk.columns[(df_bulk.columns.get_level_values("condition") == condition) & (df_bulk.columns.get_level_values(level) == ko)]

                    if b.shape[0] == 0 or s.shape[0] == 0:
                        continue

                    single_cell = df[s]
                    sm = single_cell.median(axis=1).ix[gene_filter]
                    bulk = df_bulk[b]
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

                    # axis[i][j].scatter(bm, sm, alpha=0.75, s=4)

                    # axis[i][j].set_title("{} {}".format(condition, ko))
                    # axis[i][j].text(bm.max(), sm.max(), "Pearson: {0:0.3f}\n Spearman: {1:0.3f}\nRMSE: {2:0.3f}".format(*[x[0] for x in res]))
                    # axis[i][j].set_xlabel("Bulk 3' RNA-seq")
                    # axis[i][j].set_ylabel("Mean of single-cells")
            # sns.despine(fig)
            # fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.scatter.{}.png".format(level)), dpi=300, bbox_inches="tight")

            comparisons.to_csv(os.path.join("results", "bulk", "bulk_single-cell_comparison.metrics.csv"))

            comparisons2 = comparisons.copy()
            comparisons2.loc[comparisons2['metric'] == "n_cells", "value"] = np.log10(comparisons2.loc[comparisons2['metric'] == "n_cells", "value"])

            # heatmap with metrics
            fig, axis = plt.subplots(
                len(comparisons2["metric"].drop_duplicates()),
                len(df.columns.get_level_values('condition').drop_duplicates()),
            )
            for i, condition in enumerate(df.columns.get_level_values('condition').drop_duplicates()):
                for j, metric in enumerate(comparisons2["metric"].drop_duplicates()):
                    pivot = pd.pivot_table(comparisons2[
                        (comparisons2["gene_group"] == gene_group) &
                        (comparisons2["level"] == level) &
                        (comparisons2["condition"] == condition) &
                        (comparisons2["metric"] == metric)], values="value", index="metric", columns="KO")
                    sns.heatmap(pivot, ax=axis[j][i], cmap='inferno')  # , annot=True, fmt=".2f", square=True)

                    if j != 4:
                        axis[j][i].set_xticklabels(axis[j][i].get_xticklabels(), visible=False)
                        axis[j][i].set_xlabel(None, visible=False)
                    else:
                        axis[j][i].set_xticklabels(axis[j][i].get_xticklabels(), rotation=90)
            fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.{}.metric_heatmap.{}.svg".format(gene_group, level)), dpi=300, bbox_inches="tight")

            # scatter of n_cells vs metrics
            pivot = pd.pivot_table(comparisons2[
                (comparisons2["gene_group"] == gene_group) &
                (comparisons2["level"] == level)], values="value", index=["condition", "KO"], columns="metric")
            g = sns.pairplot(data=pivot.reset_index().set_index("KO"), hue="condition")
            g.fig.savefig(os.path.join("results", "bulk", "bulk_single-cell_comparison.{}.metric_scatter.{}.svg".format(gene_group, level)), dpi=300, bbox_inches="tight")

            comparisons.to_csv(os.path.join("results", "bulk", "bulk_single-cell_comparison.metrics.csv"))

    # 2b).
    # CD69, CD82, etc...
    # At gene level or at grna level
    for level in ["gene", "grna"]:
        print(level)
        markers = ["CD69", "CD82", "PDCD1", "CD38", "BCL7A", "CDC20", "TUBB", "ADA", "TUBA1B"]

        sm = df.T.groupby(level=['condition', level]).mean()[markers]
        sm["data_type"] = "single_cell"
        bm = df_bulk.T.groupby(level=['condition', level]).mean()[markers]
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
            len(df.columns.get_level_values('condition').drop_duplicates()),
            len(markers),
            figsize=(
                4 * len(markers),
                4 * len(df.columns.get_level_values('condition').drop_duplicates()))
        )
        for i, condition in enumerate(df.columns.get_level_values('condition').drop_duplicates()):
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
    for i, condition in enumerate(df_bulk.columns.levels[1]):
        for j, ko in enumerate(df_bulk.columns.levels[2]):
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

    for i, condition in enumerate(df_bulk.columns.levels[1]):
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
    for i, condition in enumerate(df_bulk.columns.levels[1]):
        for j, ko in enumerate(df_bulk.columns.levels[2]):
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
    for i, condition in enumerate(df_bulk.columns.levels[1]):
        for j, ko in enumerate(df_bulk.columns.levels[2]):
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

    for i, condition in enumerate(df_bulk.columns.levels[1]):
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


def flow_analysis(df, flow_df):
    """
    """
    from FlowCytometryTools import FCMeasurement
    from FlowCytometryTools import ThresholdGate, PolyGate
    from scipy.stats import pearsonr, spearmanr

    #

    # Plot pregated mean of positive cells vs marker gene expression per gRNA
    marker_gene_mapping = {
        "CD38": "CD38", "CD69": "CD69", "PD1": "PDCD1", "pRelA": "RELA", "CD154": "CD40LG", "CD25": "IL2RA", "CD82": "CD82"
    }
    df_t = df.T

    # Each marker separately, both conditions together
    condition_dict = dict(zip(flow_df['condition'].drop_duplicates().tolist()[::-1], sns.color_palette("colorblind") * 10))

    for marker in marker_gene_mapping.keys():
        for level in ["gene", "grna"]:
            for metric in ["%+cells", "gMFI"]:
                fig, axis = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(1 * 4, 1 * 4))

                facs_means = dict()
                scRNA_means = dict()
                for i, condition in enumerate(flow_df['condition'].drop_duplicates()):
                    for level_value in flow_df[flow_df['condition'] == condition][level].drop_duplicates():
                        if metric == "%+cells":
                            facs_mean = flow_df[(flow_df['condition'] == condition) & (flow_df[level] == level_value)]["%{}+".format(marker)].squeeze().mean()
                        else:
                            facs_mean = flow_df[(flow_df['condition'] == condition) & (flow_df[level] == level_value)]["{} gMFI".format(marker)].squeeze().mean()
                        if (level == "grna") and (~level_value.startswith("CTRL")):
                            level_value = "Tcrlibrary_" + level_value
                        scRNA_mean = df_t.ix[
                            df_t.index[
                                (df_t.index.get_level_values("condition") == condition) &
                                (df_t.index.get_level_values(level) == level_value)
                            ]][marker_gene_mapping[marker]].mean()
                        if pd.isnull(facs_mean) or pd.isnull(scRNA_mean):
                            continue

                        axis.scatter(scRNA_mean, facs_mean, s=10, alpha=0.7, color=condition_dict[condition])
                        axis.text(scRNA_mean, facs_mean, level_value)
                        facs_means[(condition, level_value)] = facs_mean
                        scRNA_means[(condition, level_value)] = scRNA_mean
                    axis.set_title("{}".format(marker))
                    axis.set_xlabel("{} expression log2(1 + TPM)".format(marker))
                    axis.set_ylabel("% cells expressing {} (FACS)".format(marker))
                if (len(scRNA_means) == 0) or (len(facs_means) == 0):
                    continue
                p, _ = pearsonr(scRNA_means.values(), facs_means.values())
                s, _ = spearmanr(scRNA_means.values(), facs_means.values())
                axis.text(max(scRNA_means.values()), min(facs_means.values()), "Pearson = {0:0.3f}\nSpearman = {1:0.3f}".format(p, s))
                sns.despine(fig)
                fig.savefig(os.path.join("results", "flow", "flow.all_samples.mean_flow_mean_scRNA.{}.{}_only.both_conditions.{}.svg".format(level, marker, metric)), bbox_inches="tight")

    # parse fcs files,
    # gate on single, live cells
    # extract cells passing crteria
    fig0, axis0 = plt.subplots(8, 12, sharex=True, sharey=True, figsize=(12 * 2, 8 * 2))
    fig1, axis1 = plt.subplots(8, 12, sharex=True, sharey=True, figsize=(12 * 2, 8 * 2))
    fig2, axis2 = plt.subplots(8, 12, sharex=True, sharey=True, figsize=(12 * 2, 8 * 2))
    fcs_cells = dict()
    for i, fcs_file in enumerate(flow_df['fcs_file']):
        if (flow_df.loc[flow_df['fcs_file'] == fcs_file, "failed"]).squeeze() is True:
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
        condition = (flow_df.loc[flow_df['fcs_file'] == fcs_file, "condition"]).squeeze()
        grna = (flow_df.loc[flow_df['fcs_file'] == fcs_file, "grna"]).squeeze()
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


def main():

    prj = Project(os.path.join("metadata", "config.yaml"))
    prj.add_sample_sheet()
    prj.paths.results_dir = results_dir = os.path.join("results")

    # get guide annotation
    guide_annotation = os.path.join("metadata", "guide_annotation.csv")
    guide_annotation = pd.read_csv(guide_annotation)

    # get expression
    for experiment in prj.sheet.df['experiment'].dropna().drop_duplicates():
        for n_genes in [500]:
            print(experiment, n_genes)

            # Read in digital expression and add pd.MultiIndex with metadata
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

            # Normalize data
            matrix_norm = normalize(exp_assigned, experiment=experiment, kind="total")
            matrix_norm = pd.read_hdf(os.path.join(results_dir, "digital_expression.500genes.{}.log2_tpm.hdf5.gz".format(experiment)), "log2_tpm", compression="gzip")

            # Filtering/Cleanup
            # remove lowly expressed genes
            df = matrix_norm[matrix_norm.sum(1) > 10]
            # remove gRNAs
            df = df[~df.index.str.contains("library|CTRL")]
            # remove ribosomal, mitochondrial genes
            df = df[~((df.index.str.contains("^RP")) | (df.index.str.contains("^MRP")) | (df.index.str.contains("^MT-")))]
            # remove gRNAs targeting essential genes
            df = df[df.columns[~df.columns.get_level_values("gene").isin(["DHODH", "MVD", "TUBB"])]]
            df.to_hdf(os.path.join(results_dir, "digital_expression.500genes.{}.log2_tpm.filtered.hdf5.gz".format(experiment)), "log2_tpm", compression="gzip")
            df = pd.read_hdf(os.path.join(results_dir, "digital_expression.500genes.{}.log2_tpm.filtered.hdf5.gz".format(experiment)), "log2_tpm", compression="gzip")

            # Load bulk RNA-seq data
            bitseq = pd.read_csv(os.path.join("results", "{}.count_matrix.gene_level.csv".format(experiment)), index_col=[0], header=range(4))
            df_bulk = np.log2(1 + bitseq.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)
            df_bulk = df_bulk.ix[df.index].dropna()  # match quantifyied genes with CROP-seq

            # Unsupervised analysis
            # apply dimentionality reduction methods/clustering
            # and discover biological axis related with stimulation
            unsupervised(df, experiment=experiment)
            # recover genes most associated with it
            diff = differential_genes(df, method="pca", experiment=experiment)
            de_genes = diff[abs(diff) > np.percentile(abs(diff), 99)].index.tolist()
            # examine their enrichments
            enrich_signature(experiment=experiment, n_genes=n_genes, method="pca")

            # investigate intra-gene, inter-grna variability transcriptome-wide or in signature only
            intra_variability(df, df_bulk, de_genes, experiment=experiment)

            # Test significance of perturbations transcriptome-wide or in signature only
            significant_perturbation(df, df_bulk, diff, experiment=experiment)

            # visualize signature and make signature position assignments per cell/gRNA/gene
            stimulation_signature(df, df_bulk, de_genes, experiment=experiment)

            # Compare with bulk RNA-seq
            # read in diff genes from bulk
            degs = pd.read_csv(os.path.join(results_dir, "bulk", "deseq_bitseq_CTRL_samples_stimulation_signature.diff.csv"), index_col=0)
            de_genes_bulk = degs[degs["padj"] < 0.05].index.tolist()
            # first inspect bulk data
            inspect_bulk(df, df_bulk, de_genes, de_genes_bulk)
            # now compare
            compare_bulk(df, experiment="CROP-seq_Jurkat_TCR")

            # Compare with FACS measurements
            # read in FACS data
            flow_df = pd.read_csv(os.path.join("metadata", "flow_analysis.csv"))
            flow_df = flow_df[flow_df["failed"] != True]
            flow_analysis(df, flow_df)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
