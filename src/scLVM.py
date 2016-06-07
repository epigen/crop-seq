
import os
import pandas as pd
import numpy as np
import h5py
import scipy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
from scLVM import scLVM
from scLVM.utils.barplot import var_plot
from scLVM.utils.misc import PCA


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


scLVM_BASE = 'scLVM'
root_dir = "/home/afr/workspace/crop-seq/"
output_dir = os.path.join(root_dir, "quantification")

# Load R annotations
f = h5py.File('cellcycle_noise.h5f', "r")
tech_noise = f['tech_noise'][:]  # technical noise
cc_genes = f['cc_genes'][:]  # cell cycle genes
ht_genes = f['gene_names_het'][:]  # heterogeneous genes

#  Load data
exp = pd.read_csv(os.path.join(root_dir, "data/runs/BSF_0222_HEK293T_Cas9/digital_expression.100genes.assigned.tsv"), sep="\t", index_col=0)
cells_per_gene = exp.apply(lambda x: (x > 1).sum(), axis=1)
genes_per_cell = exp.apply(lambda x: (x > 1).sum(), axis=0)

# Gene/cells thresholds
n_cells = 10
n_genes = 500
matrix = exp[cells_per_gene > n_cells][exp.columns[genes_per_cell >= n_genes]]

# Log2
matrix = np.log2(1 + matrix)

# Prepare indices
geneID = matrix.index.tolist()
Y = matrix.values
genes_het_bool = matrix.index.isin(ht_genes)  # index of heterogeneous genes
idx_cell_cycle = matrix.index.isin(cc_genes)  # index of cell cycle genes
tech_noise = tech_noise[(cells_per_gene > n_cells).values]

# subset gene matrixression matrix
# filter cell cycle genes
Ycc = matrix[~matrix.index.isin(cc_genes)]
# remove genes with zero counts
Ycc = Ycc[(Ycc.mean(axis=1) ** 2) > 0]

# visualize
fig, axis = plt.subplots(1, figsize=(8, 8))
axis.imshow(Ycc.T, cmap=cm.RdBu, vmin=0, vmax=+3, interpolation='None')
axis.set_xticks([])
axis.set_yticks([])
axis.set_xlabel('genes')
axis.set_ylabel('cells')

locations = matrix.reset_index()[~matrix.index.isin(Ycc.index)].index.tolist()  # numerical index of rows not in Ycc
k = 80  # number of latent factors
out_dir = os.path.join(scLVM_BASE + 'cache')  # folder where results are cached
file_name = 'Kcc.hdf5'  # name of the cache file
recalc = True  # recalculate X and Kconf
use_ard = True  # use automatic relevance detection
sclvm = scLVM(Y.T)
# Fit model with 80 factors
X_ARD, Kcc_ARD, varGpltVM_ARD = sclvm.fitGPLVM(
    idx=locations, k=k, out_dir=out_dir, file_name=file_name, recalc=recalc, use_ard=use_ard)

# Plot variance contributions from ARD
fig, axis = plt.subplots(1)
axis.scatter(scipy.arange(k) + 1, varGpltVM_ARD['X_ARD'])
axis.set_title('Variance explained by latent factors')
axis.set_xlim([0, k + 1])
axis.set_xlabel('# Factor')
axis.set_ylabel('Variance explained')

# Heatmap
sns.clustermap(Kcc_ARD)

# Fit with one latent factor only
sclvm = scLVM(Y.T)
X, Kcc, varGpltVM = sclvm.fitGPLVM(
    idx=locations, k=1, out_dir=out_dir, file_name=file_name, recalc=True, use_ard=False)

# Heatmap
sns.clustermap(Kcc)

# Variance decomposition and cell cycle correction
# considers only heterogeneous genes
Yhet = matrix[matrix.index.isin(ht_genes)]
geneID_het = Yhet.index.tolist()

# optionally: restrict range for the analysis
i0 = 0  # gene from which the analysis starts
i1 = len(Y)  # 1000  # gene at which the analysis ends

# construct sclvm object
# sclvm = scLVM(Yhet.T.values, geneID=geneID, tech_noise=tech_noise)
sclvm = scLVM(Y.T, geneID=geneID, tech_noise=tech_noise)
# fit the model
sclvm.varianceDecomposition(K=Kcc)  # ,i0=i0,i1=i1

normalize = True  # variance components are normalizaed to sum up to one

# get variance components
var, var_info = sclvm.getVarianceComponents(normalize=normalize)
var_filtered = var[var_info['conv']]  # filter out genes for which vd has not converged

# get corrected expression levels
Ycorr = sclvm.getCorrectedExpression()
pd.DataFrame(Ycorr.T, index=matrix.index, columns=matrix.columns).to_csv(
    os.path.join(output_dir, "corrected_values.csv"), index=True)
Ycorr.shape

# Compare before/after
guides = pd.Series(matrix.columns.str.split("_")).apply(lambda x: x[1])
color_map = dict(zip(list(set(guides)), sns.color_palette("colorblind")))

sns.clustermap(Y[i0:i1, :], xticklabels=False, yticklabels=False, col_colors=[color_map[g] for g in guides])
plt.savefig(os.path.join(output_dir, "expression.png"), bbox_inches="tight", dpi=300)
plt.close()
sns.clustermap(Ycorr[i0:i1, :].T, xticklabels=False, yticklabels=False, col_colors=[color_map[g] for g in guides])
plt.savefig(os.path.join(output_dir, "expression.corrected.png"), bbox_inches="tight", dpi=300)
plt.close()


corrected = pd.DataFrame(Ycorr.T, index=matrix.index, columns=matrix.columns)
genes = ['MBD1_gene', 'DNMT3B_gene', 'filler_gene', 'TET2', 'TET2_gene']
sns.clustermap(corrected.ix[genes].dropna(), yticklabels=True, xticklabels=False, col_colors=[color_map[g] for g in guides])
plt.savefig(os.path.join(output_dir, "expression.corrected.guide_chroms.png"), bbox_inches="tight", dpi=300)
plt.close()

# maybe replace negatives with 0
# Ycorr[Ycorr < 0] = 0

# calculate average variance components across all genes and visualize
var_mean = var_filtered.mean(0)
colors = ['Green', 'MediumBlue', 'Gray']
pp = plt.pie(var_mean, labels=var_info['col_header'], autopct='%1.1f%%', colors=colors, shadow=True, startangle=0)

# visualize this stratifying for different levels of technical noise
H2 = 1 - var_filtered[:, 2]
var_comp_fileds = scipy.array([
    [0, 'cell cycle', 'Peru'],
    [1, 'biol. var', 'DarkMagenta'],
    [2, 'tech. var', '#92c5de']], dtype=object)
var_plot(var_filtered, H2, var_comp_fileds, normalize=True, figsize=[5, 4])


# Going forward

i0 = 0     # gene from which the analysis starts
i1 = 10    # gene to which the analysis ends

# fit lmm without correction
pv0, beta0, info0 = sclvm.fitLMM(K=None, i0=i0, i1=i1, verbose=False)
# fit lmm with correction
pv1, beta1, info1 = sclvm.fitLMM(K=Kcc, i0=i0, i1=i1, verbose=False)

#

fig, axis = plt.subplots(2)
axis.flatten()
axis[0].set_title('Without Correction')
p = axis[0].imshow(beta0[:, i0:i1], cmap=cm.RdBu, vmin=-0.6, vmax=+1, interpolation='None')
plt.colorbar()
axis[0].set_xticks([])
axis[0].set_yticks([])
axis[0].set_xlabel('gene'), axis[0].set_ylabel('gene')
axis[1].set_title('With Correction')
p = axis[1].imshow(beta1[:, i0:i1], cmap=cm.RdBu, vmin=-0.6, vmax=+1, interpolation='None')
plt.colorbar()
axis[1].set_xticks([])
axis[1].set_yticks([])
axis[1].set_xlabel('gene'), axis[1].set_ylabel('gene')


# Get "differential"
diff = np.log2(
    (1 + corrected[corrected.columns[guides == "DNMT3B"]].median(1)) /
    (1 + corrected[corrected.columns[guides == "filler"]].median(1))).sort_values()

sns.clustermap(
    # corrected.ix[diff.head(100).index.tolist() + diff.tail(100).index.tolist()],
    corrected.ix[diff.tail(100).index.tolist()][corrected.columns[corrected.columns.str.contains("MBD1|filler")]],
    col_colors=np.array([color_map[g] for g in guides])[corrected.columns.str.contains("MBD1|filler")],
    standard_scale=1)


fig, axis = plt.subplots(2, 2, figsize=(12, 8))
axis = axis.flatten()
for i, gene in enumerate(["MBD1", "TET2", "DNMT3B"]):
    axis[i].scatter(
        corrected[corrected.columns[guides == gene]].median(1),
        corrected[corrected.columns[guides == "filler"]].median(1)
    )
    axis[i].set_title(gene)
sns.despine(fig)
fig.savefig(os.path.join(output_dir, "expression.corrected.comparison_guides.png"), bbox_inches="tight", dpi=300)


corrected[
    ((corrected[corrected.columns[guides == "MBD1"]].median(1)) > 1.5) &
    ((corrected[corrected.columns[guides == "filler"]].median(1)) < 1)
]

import GPy
m = Ycorr.copy()
Ycorr = Ycorr.T.values
# Model optimization
Ystd = Ycorr - Ycorr.mean(0)
Ystd /= Ystd.std(0)
input_dim = 2  # How many latent dimensions to use
kern = GPy.kern.RBF(input_dim, ARD=True)  # ARD kernel
m = GPy.models.BayesianGPLVM(Ystd, input_dim=input_dim, kernel=kern, num_inducing=40)
m.optimize('scg', messages=0, max_iters=2000)


m.kern.plot_ARD()

# color mapping
color_dict = dict(zip(
    set(pd.Series(Ycorr.index.str.split("_")).apply(lambda x: x[1])),
    sns.color_palette("colorblind")
))
colors = [color_dict[q] for q in pd.Series(Ycorr.index.str.split("_")).apply(lambda x: x[1])]
fig, axis = plt.subplots(2)

# linear PCA
from scLVM.utils.misc import PCA
[S, W] = PCA(Ystd, 2)
axis[0].scatter(S[:, 0], S[:, 1], 40, colors)
axis[0].set_xlabel('PC1')
axis[0].set_ylabel('PC2')
# Non-linear PCA
axis[1].scatter(m.X[:, 0].mean, m.X[:, 1].mean, 40, colors)
axis[1].set_xlabel('PC1')
axis[1].set_ylabel('PC2')
fig.savefig(os.path.join(output_dir, "clustering.pca.corrected.png"), bbox_inches="tight", dpi=300)


# Inspect
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS, LocallyLinearEmbedding, SpectralEmbedding, Isomap

fig, axis = plt.subplots(3, 2, figsize=(12, 8))
axis = axis.flatten()
colors = [color_dict[q] for q in pd.Series(Ycorr.index.str.split("_")).apply(lambda x: x[1])]
for i, method in enumerate([PCA, MDS, TSNE, LocallyLinearEmbedding, SpectralEmbedding, Isomap]):
    fitted = method().fit_transform(Ycorr)
    axis[i].scatter(fitted[:, 0], fitted[:, 1], color=colors)
fig.savefig(os.path.join(output_dir, "clustering.corrected.png"), bbox_inches="tight", dpi=300)
