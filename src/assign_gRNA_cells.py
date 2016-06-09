
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


def get_reads_in_construct(bam, guide_annotation):
    def overlap_1d(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    # sequence templates
    u6 = (
        "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATAC"
        "GATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGAC"
        "TGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGA"
        "AAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTAT"
        "GTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAA"
        "GTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC")
    g = "G"
    "CAGGATTGGGGGCGAGTCGG"  # DNMT3B as example
    rest = (
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTA"
        "GTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGT"
        "GCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACAC"
        "TGCTTTTTGCTTGTACTGGGTCTCTCTGGTTAGACCA"
        "GATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAAC"
        "CCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGC"
        "TTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGG"
        "TAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGG"
        "AAAATCTCTAGCAGTACGTATAGTAGTTCATGTCATC"
        "TTATTATTCAGTATTTATAACTTGCAAAGAAATGAAT"
        "ATCAGAGAGTGAGAGGAACTTGTTTATTGCAGCTTAT"
        "AATGGTTACAAATAAAGCAATAGCATCACAAATTTCA"
        "CAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGG"
        "TTTGTCCAAACTCATCAATGTATCTTATCATGTCTG")

    bam_handle = pysam.AlignmentFile(bam)

    reads = pd.DataFrame()
    # for each "chromosome" (each guideRNA)
    for chrom in guide_annotation["target"].unique():
        alns = bam_handle.fetch(reference=chrom + "_chrom")

        # get position of alignment
        guide_seq = guide_annotation[guide_annotation["target"] == chrom]['sequence'].squeeze()
        chrom_size = len(u6 + g + guide_seq + rest)
        guide_start_pos = len(u6) + 1
        guide_end_pos = chrom_size - len(rest)

        # for each read
        for aln in alns:
            # skip reads
            if (
                aln.is_qcfail or  # failed quality (never happens, but for the future)
                np.mean(aln.query_alignment_qualities) < 10 or  # low mapping Q (never happens, but for the future)
                "--" in aln.get_reference_sequence()  # reads with two+ gaps
            ):
                continue

            # get cell index
            cell = dict(aln.get_tags())['XC']

            # get molecule index
            molecule = dict(aln.get_tags())['XM']

            # determine if inside gRNA
            if (aln.reference_end >= guide_start_pos) and (aln.reference_start < guide_end_pos):
                inside = True
            else:
                inside = False
            # determine distance to end of gRNA sequence
            distance = aln.reference_start - guide_end_pos

            # determine numbner of overlaping bases
            overlap = overlap_1d(aln.reference_start, aln.reference_end, guide_start_pos, guide_end_pos)

            # make sure strand is correct
            if aln.is_reverse:
                strand_agreeement = False
            else:
                strand_agreeement = True
            # get alignement quality
            # aln.mapapping_quality
            mapping_quality = np.mean(aln.query_alignment_qualities)

            reads = reads.append(pd.Series([chrom, cell, molecule, distance, overlap, inside, mapping_quality, strand_agreeement]), ignore_index=True)
    reads.columns = ["chrom", "cell", "molecule", "distance", "overlap", "inside", "mapping_quality", "strand_agreeement"]
    return reads


def plot_reads_in_constructs(reads):
    from collections import Counter
    # Inspect
    fig, axis = plt.subplots(2)
    # number of barcode reads per cell
    sns.distplot(np.log2(reads.groupby(["cell"]).apply(len)), ax=axis[0], kde=False)
    # number of unique barcode reads per cell
    sns.distplot(np.log2(reads.groupby(["cell"])['molecule'].apply(set).apply(len)), ax=axis[1], kde=False)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.svg"), bbox_inches="tight")

    # index
    u = reads.ix[reads[["cell", 'molecule']].drop_duplicates().index]

    # number of unique barcode reads saying inside per cell
    fig, axis = plt.subplots(1)
    sns.distplot(np.log2(u.groupby(["cell"])['molecule'].apply(len)), ax=axis, kde=False)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.inside.svg"), bbox_inches="tight")

    # concordance between reads in same cell
    inside_ratio = u.groupby(["cell"]).apply(lambda x: x['inside'] / len(x))
    fig, axis = plt.subplots(1)
    sns.distplot(inside_ratio, kde=False)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.concordance.svg"), bbox_inches="tight")

    # how many cells are fully out?
    inside_ratio[inside_ratio == 0.0].shape[0]
    # inside?
    inside = inside_ratio[1.0 == inside_ratio]
    inside.shape[0]
    # in-between?
    between = inside_ratio[(1.0 > inside_ratio) & (inside_ratio > 0.0)]
    between.shape[0]

    # distribution of guides for cells inside
    fig, axis = plt.subplots(2)
    c = Counter(u[u['cell'].isin(inside.reset_index()['cell'])][['chrom', 'cell']].drop_duplicates()['chrom'])
    sns.barplot(c.keys(), c.values(), ax=axis[0])
    c = Counter(u[u['cell'].isin(between.reset_index()['cell'])][['chrom', 'cell']].drop_duplicates()['chrom'])
    sns.barplot(c.keys(), c.values(), ax=axis[1])
    axis[0].set_title("Reads inside")
    axis[0].set_title("Reads between")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.guides_inside.svg"), bbox_inches="tight")

    # distribution of reads regarding constructs (for each guide)
    g = sns.FacetGrid(u, col="chrom", sharex=False, sharey=False)
    g.map(sns.distplot, 'distance')
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.distance.svg"), bbox_inches="tight")

    # Inspect weird MBD1 distribution pattern
    fig, axis = plt.subplots(2)
    u[u['chrom'] == "MBD1"].groupby("inside")['distance'].apply(sns.distplot, ax=axis[0])
    u[u['chrom'] == "MBD1"].groupby("inside")['overlap'].apply(sns.distplot, ax=axis[1])
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.guides_inside.MBD1_distance_overlap_in_out.svg"), bbox_inches="tight")
    # Inspect weird DNMT3B distribution pattern
    fig, axis = plt.subplots(2)
    u[u['chrom'] == "DNMT3B"].groupby("inside")['distance'].apply(sns.distplot, ax=axis[0])
    u[u['chrom'] == "DNMT3B"].groupby("inside")['overlap'].apply(sns.distplot, ax=axis[1])
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.guides_inside.DNMT3B_distance_overlap_in_out.svg"), bbox_inches="tight")

    # For all guides plot distribution inside and out
    g = sns.FacetGrid(u, col="chrom", row="inside", sharex=False, sharey=False)
    g.map(sns.distplot, 'distance')
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.distance.in_out.svg"), bbox_inches="tight")

    g = sns.FacetGrid(u, col="chrom", row="inside", sharex=False, sharey=False)
    g.map(sns.distplot, 'overlap')
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.overlap.in_out.svg"), bbox_inches="tight")


def make_assignment(reads, cells):
    def score(cell):
        return cell.groupby("chrom")['overlap'].apply(lambda x: x.sum() / len(x))

    # Assign
    # unique reads per cell
    u = reads.ix[reads[["cell", 'molecule']].drop_duplicates().index]
    # filter out reads in wrong strand
    u = u[u['strand_agreeement'] == 1]
    # get unique reads per cell
    u = reads.ix[reads[["cell", 'molecule']].drop_duplicates().index]

    scores = u.groupby(["cell"]).apply(score)
    scores = scores.reset_index().pivot_table("overlap", "cell", "chrom").fillna(0)

    # assign (get max)
    scores["assignment"] = scores.apply(np.argmax, axis=1)
    scores["score"] = scores[list(set(reads['chrom']))].apply(max, axis=1)
    # give nan to cells with no overlap
    scores.loc[scores['score'] == 0, 'assignment'] = pd.np.nan
    scores.loc[scores['score'] == 0, 'score'] = pd.np.nan

    # concordance between reads in same cell
    concordance_ratio = scores.drop(["assignment", "score"], axis=1).apply(lambda x: x / sum(x), axis=1).replace({0.0: pd.np.nan})

    g = sns.FacetGrid(pd.melt(concordance_ratio, value_name="score_concordance"), col="chrom", sharex=True, sharey=False)
    g.map(sns.distplot, 'score_concordance', kde=False)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.score_concordance_ratio.svg"), bbox_inches="tight")

    # join data from all cells for completeness
    cells = pd.Series(dict(zip(cells, [pd.np.nan] * len(cells)))).reset_index()
    cells.columns = ["cell", "assignment"]
    cells["score"] = pd.np.nan
    assignments = scores.reset_index()[["cell", "assignment", "score"]]
    complete_assignment = pd.merge(assignments, cells[~cells['cell'].isin(assignments['cell'])], how="outer")

    # plot assignment stats
    assignments.dropna().shape[0]
    c = assignments['assignment'].value_counts()  # .replace(pd.np.nan, "None")
    fig, axis = plt.subplots(1)
    sns.barplot(c.index, c.values, ax=axis)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.identified.svg"), bbox_inches="tight")

    # plot assignment stats
    c = complete_assignment['assignment'].value_counts()  # .replace(pd.np.nan, "None")
    fig, axis = plt.subplots(1)
    sns.barplot(c.index, c.values, ax=axis)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.identified.abstotal.svg"), bbox_inches="tight")

    fig, axis = plt.subplots(2, 2, figsize=(8, 8), sharex=False, sharey=False)
    axis = axis.flatten()
    for i, gene in enumerate(set(scores['assignment'].dropna())):
        sns.distplot(scores[scores['assignment'] == gene]['score'], kde=False, ax=axis[i])
        axis[i].set_title(gene)
        axis[i].set_xlim(0, max(scores[scores['assignment'] == gene]['score']))
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.scores.distibution.svg"), bbox_inches="tight")

    return scores, complete_assignment


def annotate_expression_with_guide(expression, complete_assignment, suffix="assigned", save=True):
    import re
    output_name = re.sub("tsv", "%s.tsv" % suffix, expression)

    exp = pd.read_csv(os.path.join(root_dir, expression), sep="\t", index_col=0)
    cm = complete_assignment.set_index('cell')['assignment'].to_dict()
    exp.columns = ["_".join([cell, str(cm[cell])]) for cell in exp.columns]
    # save
    if save:
        exp.to_csv(output_name, sep="\t", index=True)
    return exp


def big_heatmap(x, complete_assignment):
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
    selected_cells = complete_assignment.set_index("cell").ix[pd.Series(x_plot.columns.str.split("_")).apply(lambda x: x[0])]
    integer_map = dict([(val, i) for i, val in enumerate(set(pd.Series(exp.columns.str.split("_")).apply(lambda x: x[1])))])
    selected_cells["assignment_color"] = selected_cells.apply(lambda x: integer_map[str(x['assignment'])], axis=1)

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
    im2 = ax2.imshow(selected_cells[["assignment_color"]].values.T, interpolation='nearest', aspect='auto', cmap="viridis")
    ax2.set_xticks([])
    ax2.set_yticks([])

    # plot gRNA assignment score
    ax3 = fig.add_axes([0.1, 0.925, 0.8, 0.015])
    im3 = ax3.imshow(np.log2(selected_cells[["score"]].values.T), interpolation='nearest', aspect='auto')
    ax3.set_xticks([])
    ax3.set_yticks([])

    fig.colorbar(
        im1, ax=ax1, label='log2(expression)',
        aspect=4, shrink=0.2, orientation="horizontal", fraction=0.05, pad=0.01, anchor=(0.9, 1.0))
    # fig.colorbar(
    #     im2, ax=ax1, label='gRNA',
    #     aspect=4, shrink=0.2, orientation="horizontal", fraction=0.05, pad=0.01, anchor=(0.5, 1.0))
    # fig.colorbar(
    #     im3, ax=ax1, label='assignment score',
    #     aspect=4, shrink=0.2, orientation="horizontal", fraction=0.05, pad=0.01, anchor=(0.1, 1.0))

    return fig


# root_dir = "/scratch/users/arendeiro/dropseqtests/BSF_0222_HEK293T_Cas9/"
root_dir = "/home/afr/workspace/crop-seq/"
output_dir = os.path.join(root_dir, "quantification_500genes")

for dir in [root_dir, output_dir]:
    if not os.path.exists(dir):
        os.makedirs(dir)

# get guide annotation
# annotation = os.path.join(root_dir, "guide_annotation.csv")
annotation = os.path.join(root_dir, "metadata/guide_annotation.csv")
guide_annotation = pd.read_csv(annotation)

# read in alignments
# bam = os.path.join(root_dir, "star_gene_exon_tagged.bam")
bam = os.path.join(root_dir, "data/runs/BSF_0222_HEK293T_Cas9/star_gene_exon_tagged.bam")
reads = get_reads_in_construct(bam, guide_annotation)
reads.to_csv(os.path.join(output_dir, "guide_cell_quantification.csv"))
reads = pd.read_csv(os.path.join(output_dir, "guide_cell_quantification.csv"), index_col=0)

# plot
plot_reads_in_constructs(reads)

# assign
bam_handle = pysam.AlignmentFile(bam)
cells = set([dict(aln.get_tags())['XC'] for aln in bam_handle])
scores, complete_assignment = make_assignment(reads, cells)
scores.to_csv(os.path.join(output_dir, "guide_cell_scores.csv"), index=True)
complete_assignment.to_csv(os.path.join(output_dir, "guide_cell_assignment.csv"), index=False)
scores = pd.read_csv(os.path.join(output_dir, "guide_cell_scores.csv"), index_col=0)
complete_assignment = pd.read_csv(os.path.join(output_dir, "guide_cell_assignment.csv"))


# Read in expression
# write assignement in name of cell in expression matrix
expression_matrix = os.path.join(root_dir, "data/runs/BSF_0222_HEK293T_Cas9/digital_expression.500genes.tsv")
exp = annotate_expression_with_guide(expression_matrix, complete_assignment, save=True)
# exp = pd.read_csv(os.path.join(root_dir, "data/runs/BSF_0222_HEK293T_Cas9/digital_expression.100genes.assigned.tsv"), sep="\t", index_col=0)

# stats
cells_per_gene = exp.apply(lambda x: (x > 1).sum(), axis=1)
genes_per_cell = exp.apply(lambda x: (x > 1).sum(), axis=0)

# plot number of genes assigned depending on the minimum number of genes required
total = list()
assigned = list()
for j in range(1, 5000, 10):
    pos = genes_per_cell > j
    total.append(pos.sum())
    original_names = pd.Series(genes_per_cell[pos].index.str.split("_")).apply(lambda x: x[0])
    assigned.append(original_names.isin(scores.dropna().index).sum())

fig, axis = plt.subplots(2, sharex=False, sharey=False, figsize=(8, 8))
axis[0].plot((range(1, 5000, 10)), (total), lw=4, color="b", label="All cells")
axis[0].plot((range(1, 5000, 10)), (assigned), lw=4, color="r", label="Cells with gRNA assigned")
axis[1].plot(np.log2(range(1, 5000, 10)), np.log2(total), lw=4, color="b", label="All cells")
axis[1].plot(np.log2(range(1, 5000, 10)), np.log2(assigned), lw=4, color="r", label="Cells with gRNA assigned")
axis[0].set_xlabel("Genes detected")
axis[1].set_xlabel("Genes detected (log2)")
axis[0].set_ylabel("Cells")
axis[1].set_ylabel("Cells (log2)")
plt.legend()
sns.despine(fig)
fig.savefig(os.path.join(output_dir, "barcodes_per_cell.identified.varying_genes_detected.svg"), bbox_inches="tight")


# select genes present in at least n cells
# select cells with at least n genes
for n_cells in [3, 10, 100]:
    for n_genes in [100, 500, 1000]:
        # n_cells = 100
        # n_genes = 100
        matrix = exp[cells_per_gene > n_cells][exp.columns[genes_per_cell >= n_genes]]

        # get number of cells with >500 genes per gRNA
        matrix.columns.str.contains("_nan").sum()
        matrix.columns.str.contains("_TET2").sum()
        matrix.columns.str.contains("_DNMT3B").sum()
        matrix.columns.str.contains("_MBD1").sum()
        matrix.columns.str.contains("_filler").sum()

        # transform
        matrix = np.log2(1 + matrix)

        # normalize by total
        matrix_norm = matrix.apply(lambda x: x / x.sum(), axis=0) * 1e4

        # expression heatmaps
        # all genes/cells
        # f = big_heatmap(np.log2(1 + exp).apply(lambda x: x / x.sum(), axis=0) * 1e4, complete_assignment)
        # f.savefig(os.path.join(output_dir, "expression_matrix.total.png"), bbox_inches="tight", dpi=300)
        # reduced
        for name, m in [("", matrix), ("norm-", matrix_norm)]:
            f = big_heatmap(m, complete_assignment)
            f.savefig(os.path.join(output_dir, "expression_matrix.%s%icells_%igenes.png" % (name, n_cells, n_genes)), bbox_inches="tight", dpi=300)

        # Inspect
        from sklearn.decomposition import PCA
        from sklearn.manifold import TSNE, MDS, LocallyLinearEmbedding, SpectralEmbedding, Isomap

        # color mapping
        color_dict = dict(zip(
            set(pd.Series(exp.columns.str.split("_")).apply(lambda x: x[1])),
            sns.color_palette("colorblind")
        ))

        fig, axis = plt.subplots(3, 2, figsize=(12, 8))
        axis = axis.flatten()
        colors = [color_dict[q] for q in pd.Series(exp.columns.str.split("_")).apply(lambda x: x[1])]
        for i, method in enumerate([PCA, MDS, TSNE, LocallyLinearEmbedding, SpectralEmbedding, Isomap]):
            fitted = method().fit_transform(matrix.T)
            axis[i].scatter(fitted[:, 0], fitted[:, 1], color=colors)
        fig.savefig(os.path.join(output_dir, "clustering.%icells_%igenes.png" % (n_cells, n_genes)), bbox_inches="tight")

        fig, axis = plt.subplots(3, 2, figsize=(12, 8))
        axis = axis.flatten()
        colors = [color_dict[q] for q in pd.Series(exp.columns.str.split("_")).apply(lambda x: x[1])]
        for i, method in enumerate([PCA, MDS, TSNE, LocallyLinearEmbedding, SpectralEmbedding, Isomap]):
            fitted = method().fit_transform(matrix_norm.T)
            axis[i].scatter(fitted[:, 0], fitted[:, 1], color=colors)
        fig.savefig(os.path.join(output_dir, "clustering.norm-%icells_%igenes.png" % (n_cells, n_genes)), bbox_inches="tight")

sns.clustermap(matrix_norm.ix[['TET2_gene', 'MBD1_gene', 'filler_gene']])

sns.clustermap(matrix_norm.ix[['TET2', 'TET2_gene', 'MBD1_gene', 'filler_gene']])
