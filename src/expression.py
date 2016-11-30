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

    # def quadratic_signature(array, matrix):
    #     import rpy2.robjects as robjects
    #     import rpy2.robjects.numpy2ri
    #     rpy2.robjects.numpy2ri.activate()

    #     import scipy
    #     import numpy as np
    #     from scipy import optimize

    #     """
    #     minimize
    #         F = x[1]^2 + x[2]^2 x[2]

    #     subject to:
    #         1) both have to be higher or equal to 0
    #          x[1] >= 0
    #          x[2] >= 0
    #         2) the sum of both is 1 or less
    #          x[1] + x[2] <= 1
    #         3) x[1] <= N <= x[2], where N is the observed single-cell value at each gene
    #          x[1] >= N
    #          x[2] <= N

    #     in matrix notation:
    #         F = (1/2)*x.T*H*x + c*x
    #     subject to:
    #         Ax <= b

    #     where:
    #         H = [[1, 0],  # Dmat in R notation
    #              [0, 1]]

    #         c = [0.5, 0.5]  # initial guesses of the composition?

    #         A = [[-1, 0], # constraint 1)
    #              [0, -1],
    #              [1,  1], # constraint 2)
    #              [-1, 0], # constraint 3)
    #              [0, -1],
    #              [-1, 0],
    #              [0, -1]]

    #         b = [0,0,1,N1,N1,N2,N2]

    #     """

    #     robjects.r('library("quadprog")')
    #     qp_solve = robjects.r('quadprog::solve.QP')

    #     n_samples = matrix.shape[1]
    #     m_genes = array.shape[0]

    #     dvec = np.array([0.5, 0.5])

    #     Dmat = np.zeros((n_samples, n_samples))  # scipy.sparse.eye(2).todense()
    #     np.fill_diagonal(Dmat, 1.)

    #     Amat = np.concatenate([
    #         np.array([[-1, 0], [0, -1]]),  # constraint 1)
    #         np.array([[1, 1]]),  # constraint 2)
    #         np.array([[-1, 0], [0, -1]] * m_genes)   # constraint 3)
    #     ])
    #     bvec = np.concatenate([
    #         np.array([0]),  # constraint 1)
    #         np.array([0]),  # constraint 1)
    #         np.array([1]),  # constraint 2)
    #         np.array([i for i in array.values for _ in range(2)])   # constraint 3)
    #     ])

    #     solution, value, unconstrained_solution, iterations, lagrangian, iact = np.array(
    #         qp_solve(Dmat, dvec, Amat, bvec=bvec))

    #     def loss(x, sign=1.):
    #         return sign * (0.5 * np.dot(x.T, np.dot(H, x)) + np.dot(c, x))

    #     def jac(x, sign=1.):
    #         return sign * (np.dot(x.T, H) + c)

    #     A = np.array([[1, 0, -1, 0], [0, 1, 0, -1]]).T
    #     H = Dmat
    #     b = array
    #     c = np.array([[1., 1., 0., 0.]]).T
    #     x0 = np.random.normal(0.5, 0.3, 2)

    #     cons = {'type': 'ineq',
    #             'fun': lambda x: b - np.dot(A, x),
    #             'jac': lambda x: -A}
    #     opt = {'disp': False}

    #     res_cons = optimize.minimize(loss, x0, jac=jac, constraints=cons,
    #                                  method='SLSQP', options=opt)
    #     res_uncons = optimize.minimize(loss, x0, jac=jac, method='SLSQP',
    #                                    options=opt)

    #     # Linear programming
    #     c = array.values[:4]
    #     A_ub = np.array([[1, 0, -1, 0], [0, 1, 0, -1]])
    #     A_eq = np.array([[1.], [1.]]).T
    #     b_ub = np.array([[1., 1., 0., 0.]]).T
    #     b_eq = np.array([[1.]]).T
    #     x0_bounds = matrix.values[:, 0]
    #     x1_bounds = matrix.values[:, 1]

    #     from scipy.optimize import linprog
    #     res = linprog(c, A_ub=A_ub, A_eq=A_eq, b_ub=b_ub, b_eq=b_eq, bounds=(x0_bounds, x1_bounds),
    #                   options={"disp": True})

    #     # linear programming step implemented
    #     from scipy.optimize import minimize

    #     def loss(x, sign=1.0):
    #         # return np.dot(x[0], s1) + np.dot(x[1], s2) - c
    #         return np.dot(s1, x) + c

    #     opt = {'disp': True}
    #     cons = [{
    #             'type': 'eq',
    #             'fun': lambda x: x[0] + x[1] - 1},
    #             {
    #             'type': 'ineq',
    #             'fun': lambda x: x[0] - 1},
    #             {
    #             'type': 'ineq',
    #             'fun': lambda x: x[1] - 1},
    #             {
    #             'type': 'ineq',
    #             'fun': lambda x: -x[0]},
    #             {
    #             'type': 'ineq',
    #             'fun': lambda x: -x[1]}]

    #     s1 = matrix[cond1].values
    #     s2 = matrix[cond2].values
    #     c = array.values
    #     x0 = np.random.normal(0.5, 0.3, 2)

    #     res_cons = optimize.minimize(loss, x0, method='SLSQP', constraints=cons, options=opt)

    #     # quadratic programming step implemented
    #     from scipy.optimize import minimize

    #     def loss2(x, sign=1.):
    #         return sign * (0.5 * np.dot(x.T, np.dot(H, x)) + np.dot(c, x))

    #     s1 = matrix[cond1].values
    #     s2 = matrix[cond2].values
    #     c = array.values
    #     x0 = np.random.normal(0.5, 0.3, 2)

    #     opt = {'disp': True}
    #     cons2 = [{'type': 'ineq',
    #              'fun': lambda x: b - np.dot(A, x)},
    #              {'type': 'eq',
    #              'fun': lambda x: 1 - sum(x)}]

    #     res_cons2 = optimize.minimize(loss2, x0, method='SLSQP', constraints=cons2, options=opt)

    # def quadratic_signature(array, matrix):
    #     import rpy2.robjects as robjects
    #     import rpy2.robjects.numpy2ri
    #     rpy2.robjects.numpy2ri.activate()

    #     import scipy
    #     import numpy as np
    #     from scipy import optimize

    #     # quadratic programming step implemented with new matrix
    #     from scipy.optimize import minimize

    #     def loss(x, sign=1.):
    #         return sign * (0.5 * np.dot(x.T, np.dot(H, x)) + np.dot(c, x)) - c0

    #     def jac(x, sign=1.):
    #         return sign * (np.dot(x.T, H) + c)

    #     cons = [
    #         {'type': 'ineq',
    #          'fun': lambda x: b - np.dot(A, x)},
    #         {'type': 'eq',
    #          'fun': lambda x: 1 - sum(x)}]
    #     opt = {'disp': False}

    #     H = np.array([[1, 0], [0, 1]])  # Q symmetric matrix
    #     A = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])  #
    #     b = np.array([1, 1, 0, 0])  #
    #     x0 = np.random.normal(0.5, 0.3, 2)

    #     fits_raw = list()
    #     for i in matrix.index:
    #         c = matrix.ix[i].values
    #         c0 = array.ix[i]
    #         res_cons = optimize.minimize(loss, x0, method='SLSQP', constraints=cons, options=opt)['x'][0]
    #         quadprog.solve_qp(H, c0, A, b)
    #         fits_raw.append(res_cons)
    #     return np.mean(fits_raw)

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

    # # Dimentionality reduction methods
    # # try several
    # methods = [PCA, LocallyLinearEmbedding, Isomap, SpectralEmbedding, TSNE]

    # for name, matrix in [("groups", df5), ("cells", df3.T)]:
    #     for method in methods:
    #         print(name, method.__name__)
    #         model = method()

    #         if name == "cells":
    #             color = get_grna_colors(matrix.T, assignment)
    #             m = matrix.T[matrix.dtypes == np.float64].T
    #         else:
    #             color = get_group_colors(matrix.T, assignment)
    #             m = matrix.T[matrix.dtypes == np.float64]

    #         fit = model.fit_transform(m)

    #         # plot
    #         if method.__name__ == "PCA":
    #             pcs = 3
    #         else:
    #             pcs = 1

    #         fig, axis = plt.subplots(2, pcs, sharex=False, sharey=False, figsize=(8, 10))
    #         if method.__name__ != "PCA":
    #             axis = [[x] for x in axis]

    #         for i, variable in enumerate(["condition", "gene"]):
    #             for pc in range(pcs):
    #                 axis[i][pc].scatter(fit[:, pc], fit[:, pc + 1], color=color[i], alpha=0.75 if name == "groups" else 0.1)
    #                 axis[i][pc].set_xticklabels([])
    #                 axis[i][pc].set_yticklabels([])
    #                 if method.__name__ == "PCA":
    #                     axis[i][pc].set_xlabel("PC %i" % pc)
    #                     axis[i][pc].set_ylabel("PC %i" % (pc + 1))
    #             # axis[i][pc].legend(
    #             #     handles=[mpatches.Patch(color=v, label=k) for k, v in color_mapping.items()],
    #             #     ncol=2 if feature == "organ" else 1,
    #             #     loc='center left',
    #             #     bbox_to_anchor=(1, 0.5))
    #         fig.savefig(os.path.join(results_dir, "differential_expression.{}.{}.{}.png".format(prefix, name, method.__name__)), bbox_inches="tight", dpi=300)

    #

    # Signature-based cell assignemnt

    # Get stimulation signature
    # 1. get mean expression of each group in signature gens
    df_grna_means = df.T.groupby(level=["condition", "grna"]).mean().T
    df_gene_means = df.T.groupby(level=["condition", "gene"]).mean().T
    x1 = df_gene_means[df_gene_means.columns[df_gene_means.columns.get_level_values("condition") == cond1]].mean(axis=1).ix[de_genes]
    x1.name = cond1
    x2 = df_gene_means[df_gene_means.columns[df_gene_means.columns.get_level_values("condition") == cond2]].mean(axis=1).ix[de_genes]
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

    # # Try quadratic programming
    # fits_raw = list()
    # for i, cell in enumerate(df.columns):
    #     if i % 100 == 0:
    #         print(i)
    #     fit = quadratic_signature(array=df.ix[de_genes][cell], matrix=pd.DataFrame(x1).join(x2))
    #     fits_raw.append(fit)

    # fits = pd.DataFrame(fits_raw, index=df.columns)
    # fits = fits.join(fits['x'].apply(pd.Series)).drop('x', axis=1).rename(columns={0: cond1, 1: cond2})

    # 5. investigate signatures
    sigs = cors.apply(lambda x: np.argmax(x), axis=1)
    sigs.name = "signature"
    sigs_p = p_values.apply(lambda x: np.argmax(x), axis=1)
    sigs_r = random_cors.apply(lambda x: np.argmax(x), axis=1)
    sigs_r.name = "signature"
    sigs_rp = random_p_values.apply(lambda x: np.argmax(x), axis=1)

    # get "uncorrelated" fraction (noise)
    res = 1 - cors.max(axis=1)
    res.name = "residual"
    res_r = 1 - sigs_r
    res_r.name = "residual"

    # reads per cell
    s = exp_assigned.sum(axis=0)
    s.name = "reads_per_cell"
    sigs = pd.merge(sigs.reset_index(), s.reset_index()).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])

    # get minimum correlation
    m = cors.min(axis=1)
    m.name = "min_corr"
    sigs = pd.merge(sigs.reset_index(), m.reset_index()).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])

    # get "smoothness" of ranked cumulative sum
    def ranked_cumsum(x, rank_threshold=0.25):
        # scale from 0 to 1
        xs = (x - x.min()) / (x.max() - x.min())
        rs = xs.sort_values(ascending=False).cumsum() / xs.sum()
        return (rs < rank_threshold).sum()

    r = cors.apply(ranked_cumsum, axis=1)

    # save
    cors.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_correlation.csv".format(experiment, n_genes)))
    p_values.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.matrix_p_values.csv".format(experiment, n_genes)))
    sigs.to_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation.csv".format(experiment, n_genes)))
    sigs = pd.read_csv(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation.csv".format(experiment, n_genes))).set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])
    # sigs['replicate'] = sigs['replicate'].astype(str)

    # all cells together vs random
    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4))
    sns.distplot(sigs, ax=axis[0])
    sns.distplot(p_values.apply(lambda x: np.argmin(x), axis=1), ax=axis[0])
    sns.distplot(random_cors.apply(lambda x: np.argmax(x), axis=1), ax=axis[1])
    sns.distplot(random_p_values.apply(lambda x: np.argmin(x), axis=1), ax=axis[1])
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.all_cells.correlation_pvalue.distributions.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

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

    # signature assignemnt vs residual
    g = sns.FacetGrid(
        data=sigs.reset_index(),
        col='condition')
    g.map(plt.scatter, "signature", "residual", alpha=0.2)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.correlation.signature_vs_residual.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # residual vs read coverage per cell
    g = sns.FacetGrid(
        data=sigs.reset_index(),
        row='condition', col='gene')
    g.map(plt.scatter, "signature", "residual", alpha=0.2)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.per_gene.correlation.signature_vs_residual.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    g = sns.FacetGrid(
        data=sigs.reset_index(),
        col='condition')
    g.map(plt.scatter, "signature", "reads_per_cell", alpha=0.2)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.correlation.signature_vs_readspercell.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    g = sns.FacetGrid(
        data=sigs.reset_index(),
        col='condition')
    g.map(plt.scatter, "residual", "reads_per_cell", alpha=0.2)
    g.fig.savefig(os.path.join(results_dir, "{}.{}genes.signature.correlation.residual_vs_readspercell.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)

    # plot distibution
    sigs_mean = sigs.groupby(level=['condition', 'gene']).median()
    sigs_mean["n_cells"] = sigs.groupby(level=['condition', 'gene']).apply(len)
    sigs_mean.to_csv(os.path.join(results_dir, "{}.{}genes.signature.group_means.annotated.csv".format(experiment, n_genes)), index=True)

    # Filter out genes with less than the 5th percentile of cells
    sigs_mean = sigs_mean[sigs_mean["n_cells"] >= 10]

    fig, axis = plt.subplots(1, figsize=(10, 8))
    sns.stripplot(x=sigs_mean['signature'], y=sigs_mean.index, orient="horiz", ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.all_cells.{}.mean_group_signature.strip.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Given the known condition of each group, what is the deviation from that?
    # plot as rank of mean
    fig, axis = plt.subplots(1, figsize=(10, 8))
    axis.scatter(sigs_mean['signature'].rank(ascending=False), sigs_mean['signature'])
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.all_cells.{}.mean_group_signature.rank.svg".format(experiment, n_genes)), bbox_inches="tight")

    # plot as violinplots (values per cell)
    p = sigs.reset_index().sort_values(['signature'])
    p = p[(p["gene"] != "library")]

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
    mean_matrix_norm = matrix_norm.ix[de_genes].T.groupby(level=['condition', 'gene']).mean()

    # cluster
    g = sns.clustermap(
        mean_matrix_norm.T,
        col_colors=get_level_colors(mean_matrix_norm.index),
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
    gene_order = (
        mean_matrix_norm[mean_matrix_norm.index.get_level_values('condition') == cond1].median() /
        mean_matrix_norm[mean_matrix_norm.index.get_level_values('condition') == cond2].median()).sort_values(ascending=False).index

    g = sns.clustermap(
        mean_matrix_norm[gene_order].T, z_score=0,
        col_colors=get_level_colors(mean_matrix_norm.index),
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
    mean_matrix_norm_sig_sorted = mean_matrix_norm.ix[sigs_mean.reset_index().sort_values(['condition', 'signature']).set_index(['condition', 'gene']).index]

    clust = sns.clustermap(
        mean_matrix_norm_sig_sorted,  # z_score=0,
        row_colors=get_level_colors(mean_matrix_norm_sig_sorted.index),
        metric='correlation',
        row_cluster=False, col_cluster=True,
        xticklabels=False, yticklabels=True,
        figsize=(8.62948158106742, 6))
    for item in clust.ax_heatmap.get_yticklabels():
        item.set_rotation(0)
    for item in clust.ax_heatmap.get_xticklabels():
        item.set_rotation(90)
    clust.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.gene_group_means.sorted_signature.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
    clust.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.gene_group_means.sorted_signature.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Strip plot of number of cells
    p = sigs_mean.ix[mean_matrix_norm_sig_sorted.index]

    fig, axis = plt.subplots(1, figsize=(8, 8))
    axis.scatter([1] * p.shape[0], range(p.shape[0]), s=p['n_cells'])
    [axis.text(1.0005, i, s=str(int(x))) for i, x in enumerate(p['n_cells'].values)]
    [axis.text(1.01, i, s=str(x)) for i, x in enumerate(p.index.get_level_values('condition'))]
    [axis.text(1.02, i, s=str(x)) for i, x in enumerate(p.index.get_level_values('gene'))]
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.cells_per_group.bubbles.svg".format(experiment, n_genes)), bbox_inches="tight")

    # # Single-cell matrix sorted in same way as above
    # matrix_norm.ix[de_genes].T.ix[mean_matrix_norm_sig_sorted.index]

    # p_all = df4.reset_index().set_index(["sti", "ass"])
    # # order gRNAs
    # p_all = p_all.ix[p.columns]
    # # order genes
    # p_all = p_all[[t.get_text() for t in clust.ax_heatmap.get_yticklabels()] + ["index"]]

    # # down
    # gs = stats[stats["Z"] < 0].sort_values("Z").index.tolist()
    # # up
    # gs += stats[stats["Z"] > 0].sort_values("Z", ascending=False).index.tolist()
    # p_all = p_all[gs + ["index"]]
    # # p_all = p_all[diffs.sort_values("fold_change").index.tolist() + ["index"]]

    # # for cmap in ["YlGn"]:  # , "YlOrBr", "GnBu", "Greys_r", "Oranges", "ocean_r"]:
    # g = sns.clustermap(
    #     p_all.drop("index", axis=1),
    #     # standard_scale=0,
    #     z_score=1,
    #     # cmap=cmap,
    #     # vmin=0.25,
    #     # col_colors=get_foldchange_colors(p_all.set_index("index").T, stats),
    #     row_colors=get_level_colors(p_all.set_index("index").T, assignment),
    #     metric='correlation',
    #     row_cluster=False, col_cluster=False,
    #     xticklabels=True, yticklabels=False,
    #     figsize=(6, 8.62948158106742))
    # for item in g.ax_heatmap.get_yticklabels():
    #     item.set_rotation(0)
    # for item in g.ax_heatmap.get_xticklabels():
    #     item.set_rotation(90)
    # # g.fig.savefig(os.path.join(results_dir, "differential_expression.{}.assigned_cells.sorted_signature.{}.png".format(prefix, cmap)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.assigned_cells.sorted_signature.png".format(experiment, n_genes)), bbox_inches="tight", dpi=300)
    # g.fig.savefig(os.path.join(results_dir, "{}.{}genes.differential_expression.assigned_cells.sorted_signature.svg".format(experiment, n_genes)), bbox_inches="tight")

    # Calculate deviation from CTRL for each stimulation
    m1 = sigs_mean[sigs_mean.index.get_level_values("condition") == cond1]['signature']
    diff1 = (m1 - m1[m1.index.get_level_values("gene") == "CTRL"].squeeze()).sort_values(ascending=False)
    m2 = sigs_mean[sigs_mean.index.get_level_values("condition") == cond2]['signature']
    diff2 = (m2 - m2[m2.index.get_level_values("gene") == "CTRL"].squeeze()).sort_values(ascending=False)

    fig, axis = plt.subplots(2, figsize=(4, 4 * 2), sharex=False, sharey=False)
    sns.barplot(diff1, diff1.index, orient="horiz", order=diff1.index, ax=axis[0])
    sns.barplot(diff2, diff2.index, orient="horiz", order=diff2.index, ax=axis[1])
    axis[0].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_xlabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_ylabel(cond1)
    axis[1].set_ylabel(cond2)
    # axis[0].set_xlim((0, -10))
    # axis[1].set_xlim((0, 10))
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.mean_group_signature_deviation.ranklog2.barplot.svg".format(experiment, n_genes)), bbox_inches="tight")

    fig, axis = plt.subplots(2, figsize=(4, 4 * 2), sharex=True, sharey=True)
    axis[0].scatter(diff1.rank(ascending=False), diff1)
    axis[1].scatter(diff2.rank(ascending=False), diff2)
    axis[0].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[1].set_ylabel("Signature deviation from control (propensity to cause stimulation)")
    axis[0].set_ylabel(cond1)
    axis[1].set_ylabel(cond2)
    # axis[0].set_xlim((0, -10))
    # axis[1].set_xlim((0, 10))
    sns.despine(fig)
    fig.savefig(os.path.join(results_dir, "{}.{}genes.signatures.mean_group_signature_deviation.ranklog2.scatter.svg".format(experiment, n_genes)), bbox_inches="tight")


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

n_genes = 500
experiment, rows = prj.sheet.df.groupby(['experiment']).groups.items()[1]


# get expression
for n_genes in [500]:
    for experiment, rows in prj.sheet.df.groupby(['experiment']):
        print(experiment)

        assignment = pd.read_csv(os.path.join(results_dir, "{}.guide_cell_assignment.all.csv".format(experiment)))

        # exp = pd.read_hdf(os.path.join(results_dir, "{}.digital_expression.{}genes.hdf5.gz".format(experiment, n_genes)), "exp_matrix", compression="gzip")
        exp_assigned = pd.read_hdf(os.path.join(results_dir, "{}.digital_expression.{}genes.only_assigned.hdf5.gz".format(experiment, n_genes)), "exp_matrix", compression="gzip")
        # exp_assigned = exp_assigned.T.reset_index()
        # exp_assigned['replicate'] = exp_assigned['replicate'].astype(str)
        # exp_assigned = exp_assigned.set_index(['condition', 'replicate', 'cell', 'grna', 'gene'])
        # exp_assigned = exp_assigned.T

        # exp_assigned.columns = exp_assigned.columns.set_levels([exp_assigned.columns.levels[:-1], exp_assigned.columns.levels[-1].astype(str)])

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
