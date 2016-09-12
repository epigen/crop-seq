#!/usr/bin/env python

"""
Usage:

cd /scratch/lab_bock/shared/projects/crop-seq/results_pipeline

for R in `ls`; do
sbatch -p mediumq --mem 20000 -c 4 -J assign_gRNA_cells.${R} \
-o /scratch/lab_bock/shared/projects/crop-seq/results_pipeline/${R}/assign_gRNA_cells.script.log \
~/assign_gRNA_cells.script.sh ${R}
done

"""


import argparse
import os
import pandas as pd
import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
    u6 = "".join([
        "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAAT",
        "TAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAAT",
        "TTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACT",
        "TGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG"])
    rest = "".join([
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCAC",
        "CGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTACTGG",
        "GTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTA",
        "AGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGT"])

    bam_handle = pysam.AlignmentFile(bam)

    reads = pd.DataFrame()
    # for each "chromosome" (each guideRNA)
    for chrom in guide_annotation["oligo_name"].unique():
        print(chrom)

        if chrom == "spCas9":
            continue

        # get position of alignment
        guide_seq = guide_annotation[guide_annotation["oligo_name"] == chrom]['sequence'].squeeze()
        chrom_size = len(u6 + guide_seq + rest)
        guide_start_pos = len(u6) + 1
        guide_end_pos = chrom_size - len(rest)

        # for each read
        for aln in bam_handle.fetch(reference=chrom + "_chrom"):
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

            # determine distance to end of gRNA sequence
            distance = aln.reference_start - guide_end_pos

            # determine numbner of overlaping bases
            overlap = overlap_1d(aln.reference_start, aln.reference_end, guide_start_pos, guide_end_pos)

            # determine if inside gRNA
            inside = True if overlap > 0 else False
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


def get_reads_in_Cas9_construct(bam, guide_annotation):
    def overlap_1d(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    # sequence templates
    u6 = "".join([
        "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAAT",
        "TAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAAT",
        "TTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACT",
        "TGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG"])
    rest = "".join([
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCAC",
        "CGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTACTGG",
        "GTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTA",
        "AGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGT"])

    chrom = "spCas9"

    bam_handle = pysam.AlignmentFile(bam)

    reads = pd.DataFrame()

    # get position of alignment
    guide_seq = guide_annotation[guide_annotation["target"] == chrom]['sequence'].squeeze()
    chrom_size = len(u6 + guide_seq + rest)
    guide_start_pos = len(u6) + 1
    guide_end_pos = chrom_size - len(rest)

    # for each read
    for aln in bam_handle.fetch(reference=chrom + "_chrom"):
        # skip reads
        if (
            aln.is_qcfail or  # failed quality (never happens, but for the future)
            np.mean(aln.query_alignment_qualities) < 10 or  # low mapping Q (never happens, but for the future)
            "--" in aln.get_reference_sequence()  # reads with two+ gaps
        ):
            continue

        # determine distance to start of Cas9 construct
        distance = guide_start_pos - aln.reference_start

        if distance < 0:
            continue

        # determine numbner of overlaping bases
        overlap = overlap_1d(aln.reference_start, aln.reference_end, guide_start_pos, guide_end_pos)

        # determine if overlaps Cas9 construct
        inside = True if overlap > 0 else False
        # make sure strand is correct
        if aln.is_reverse:
            strand_agreeement = False
        else:
            strand_agreeement = True
        # get alignement quality
        # aln.mapapping_quality
        mapping_quality = np.mean(aln.query_alignment_qualities)

        # get cell index
        cell = dict(aln.get_tags())['XC']

        # get molecule index
        molecule = dict(aln.get_tags())['XM']

        reads = reads.append(pd.Series([chrom, cell, molecule, distance, overlap, inside, mapping_quality, strand_agreeement]), ignore_index=True)
    reads.columns = ["chrom", "cell", "molecule", "distance", "overlap", "inside", "mapping_quality", "strand_agreeement"]
    return reads


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
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.svg"), bbox_inches="tight")

    # index
    u = reads.ix[reads[["chrom", "cell", 'molecule']].drop_duplicates().index]

    # number of unique barcode reads saying inside per cell
    fig, axis = plt.subplots(1)
    sns.distplot(np.log2(1 + u.groupby(["cell"])['molecule'].apply(len)), ax=axis, kde=False)
    axis.set_xlabel("Molecules (log2(1 + x))")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.inside.svg"), bbox_inches="tight")

    # concordance between reads in same cell
    inside_ratio = u.groupby(["cell"]).apply(lambda x: x['inside'] / len(x))
    fig, axis = plt.subplots(1)
    sns.distplot(inside_ratio, kde=False)
    axis.set_xlabel("Ratio molecules overlap gRNA / total")
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
    axis[0].set_title("Cells with only molecules inside")
    axis[1].set_title("Cells with only molecules between")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.guides_inside.svg"), bbox_inches="tight")

    # distribution of reads regarding constructs (for each guide)
    g = sns.FacetGrid(u, col="chrom", sharex=False, sharey=False)
    g.map(sns.distplot, 'distance', kde=False)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.distance.svg"), bbox_inches="tight")

    # Inspect weird MBD1 distribution pattern
    fig, axis = plt.subplots(2)
    u[u['chrom'] == "MBD1"].groupby("inside")['distance'].apply(sns.distplot, kde=False, ax=axis[0])
    u[u['chrom'] == "MBD1"].groupby("inside")['overlap'].apply(sns.distplot, kde=False, ax=axis[1])
    axis[0].set_xlabel("Distance to gRNA")
    axis[1].set_xlabel("Overlap with gRNA")
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.guides_inside.MBD1_distance_overlap_in_out.svg"), bbox_inches="tight")

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
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.distance.in_out.svg"), bbox_inches="tight")

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
                sns.distplot(subset['overlap'], kde=False, ax=ax)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.overlap.in_out.svg"), bbox_inches="tight")


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
    scores["score"] = scores[scores.columns.drop("assignment")].apply(max, axis=1)
    # give nan to cells with no overlap (this is because argmax pickus a draw)
    scores.loc[scores['score'] == 0, 'assignment'] = pd.np.nan
    scores.loc[scores['score'] == 0, 'score'] = pd.np.nan

    # Get only assigned cells
    scores = scores.dropna()

    # # Plot bp covered per cell
    # g = sns.FacetGrid(scores[["assignment", "score"]], col="assignment", sharex=False, sharey=False)
    # g.map(sns.distplot, 'score', kde=False)
    # g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.bp_covered.svg"), bbox_inches="tight")

    # # Convert to coverage in X times (divide by length of gRNA)
    # coverage = scores.drop(['assignment', 'score'], axis=1).apply(lambda x: x / float(len(guide_annotation[guide_annotation["target"] == x.name]["sequence"].squeeze())), axis=0)
    # coverage["maxscore"] = coverage.apply(max, axis=1)
    # coverage["assignment"] = coverage.drop("maxscore", axis=1).apply(np.argmax, axis=1)

    # cov = coverage[~(coverage["maxscore"] == 0)][["assignment", "maxscore"]]
    # cov["maxscore"] = np.log2(cov["maxscore"])
    # g = sns.FacetGrid(cov, col="assignment", sharex=False, sharey=False)
    # g.map(sns.distplot, 'maxscore', kde=False)
    # for a in g.axes.flatten():
    #     a.set_xlabel("coverage log2")
    # g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.coverage.svg"), bbox_inches="tight")

    # # concordance between reads in same cell
    # concordance_ratio = scores.drop(["assignment", "score"], axis=1).apply(lambda x: x / sum(x), axis=1).replace({0.0: pd.np.nan})
    # g = sns.FacetGrid(pd.melt(concordance_ratio, value_name="score_concordance"), col="chrom", sharex=True, sharey=False)
    # g.map(sns.distplot, 'score_concordance', kde=False)
    # g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.score_concordance_ratio.svg"), bbox_inches="tight")

    # Get assigned cells
    assignment = scores.reset_index()[["cell", "assignment", "score"]]

    # # plot assignment stats
    # c = assignment['assignment'].value_counts()  # .replace(pd.np.nan, "None")
    # fig, axis = plt.subplots(1)
    # sns.barplot(c.index, c.values, ax=axis)
    # fig.savefig(os.path.join(output_dir, "barcodes_per_cell.identified.svg"), bbox_inches="tight")

    # fig, axis = plt.subplots(2, 2, figsize=(8, 8), sharex=False, sharey=False)
    # axis = axis.flatten()
    # for i, gene in enumerate(set(scores['assignment'].dropna())):
    #     sns.distplot(scores[scores['assignment'] == gene]['score'], kde=False, ax=axis[i])
    #     axis[i].set_title(gene)
    #     axis[i].set_xlim(0, max(scores[scores['assignment'] == gene]['score']))
    # fig.savefig(os.path.join(output_dir, "barcodes_per_cell.scores.distibution.svg"), bbox_inches="tight")

    return scores, assignment


def annotate_expression_with_guide(expression, assignment, suffix="assigned", save=True):
    import re
    output_name = re.sub("tsv", "%s.tsv" % suffix, expression)

    exp = pd.read_csv(expression, sep="\t", index_col=0)
    cm = assignment.set_index('cell')['assignment'].to_dict()
    exp.columns = ["_".join([cell, str(cm[cell])]) if cell in cm.keys() else cell for cell in exp.columns]
    # save
    if save:
        exp.to_csv(output_name, sep="\t", index=True)
    return exp


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


# parser
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--name", dest="run_name", required=True)
args = parser.parse_args()

# dirs
root_dir = "/scratch/lab_bock/shared/projects/crop-seq"
sample_dir = os.path.join(root_dir, "results_pipeline", args.run_name)
output_dir = os.path.join(sample_dir, "quantification")

for dir_ in [root_dir, sample_dir, output_dir]:
    if not os.path.exists(dir_):
        os.makedirs(dir_)

# get guide annotation
annotation = os.path.join(root_dir, "metadata", "guide_annotation.csv")
guide_annotation = pd.read_csv(annotation)

print("Starting with sample %s" % args.run_name)

# read in alignments
bam = os.path.join(sample_dir, "star_gene_exon_tagged.bam")

# Get reads in construct
if not os.path.exists(os.path.join(output_dir, "guide_cell_quantification.csv")):
    reads = get_reads_in_construct(bam, guide_annotation)
    reads.to_csv(os.path.join(output_dir, "guide_cell_quantification.csv"), index=False)
else:
    reads = pd.read_csv(os.path.join(output_dir, "guide_cell_quantification.csv"))

# # reads in cas9 construct
# cas9_reads = get_reads_in_Cas9_construct(bam, guide_annotation)
# cas9_reads.to_csv(os.path.join(output_dir, "cas9_quantification.reads.csv"), index=False)

# cas9_expression = cas9_reads.groupby(['cell'])['molecule'].apply(np.unique).apply(len)
# cas9_expression.reset_index().to_csv(os.path.join(output_dir, "cas9_quantification.counts.csv"), index=False)

# plot
# plot_reads_in_constructs(reads)

# Assign cells to gRNA
bam_handle = pysam.AlignmentFile(bam)
cells = set([dict(aln.get_tags())['XC'] for aln in bam_handle])
scores, assignment = make_assignment(reads, cells)
scores.to_csv(os.path.join(output_dir, "guide_cell_scores.csv"), index=True)
assignment.to_csv(os.path.join(output_dir, "guide_cell_assignment.csv"), index=False)


# Read in expression
# write assignement in name of cell in expression matrix
for n_cells in [500, 100]:  # , 10]:
    if os.path.exists(os.path.join(sample_dir, "digital_expression.{}genes.tsv".format(n_cells))):
        expression_file = os.path.join(sample_dir, "digital_expression.{}genes.tsv".format(n_cells))
        exp = pd.read_csv(expression_file, sep="\t", index_col=0)
        exp = annotate_expression_with_guide(expression_file, assignment, save=True)
