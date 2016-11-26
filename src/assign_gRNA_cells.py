#!/usr/bin/env python

import os
import pandas as pd
import pysam
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


def get_reads_in_construct(bam, guide_annotation):
    def overlap_1d(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    bam_handle = pysam.AlignmentFile(bam)

    reads = pd.DataFrame()
    # for each "chromosome" (each guideRNA)
    for chrom in guide_annotation["oligo_name"].unique():
        print(chrom)

        if chrom == "Cas9_blast":
            continue

        # get position of alignment
        guide_seq = guide_annotation[guide_annotation["oligo_name"] == chrom]['sequence'].squeeze()
        chrom_size = len(prj.config['crop-seq']['u6'] + guide_seq + prj.config['crop-seq']['rest'])
        guide_start_pos = len(prj.config['crop-seq']['u6']) + 1
        guide_end_pos = chrom_size - len(prj.config['crop-seq']['rest'])

        # for each read
        for aln in bam_handle.fetch(reference=chrom + "_chrom"):
            # skip reads
            if (
                aln.is_qcfail or  # failed quality (never happens, but for the future)
                aln.is_secondary or
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

            reads = reads.append(pd.Series([
                chrom, cell, molecule, aln.reference_start, aln.reference_end,
                distance, overlap, inside, mapping_quality, strand_agreeement]), ignore_index=True)
    reads.columns = [
        "chrom", "cell", "molecule", "read_start", "read_end",
        "distance", "overlap", "inside", "mapping_quality", "strand_agreeement"]
    return reads


def get_reads_in_Cas9_construct(bam):
    def overlap_1d(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    chrom = "Cas9_blast"

    bam_handle = pysam.AlignmentFile(bam)

    reads = pd.DataFrame()

    sequence = "".join([
        prj.config['crop-seq']['cas9'],
        prj.config['crop-seq']['nls'],
        prj.config['crop-seq']['flag'],
        prj.config['crop-seq']['p2a'],
        prj.config['crop-seq']['blast'],
        prj.config['crop-seq']['space'],
        prj.config['crop-seq']['virus_ltr']])

    # get position of alignment
    chrom_size = len(sequence)
    guide_start_pos = 0
    guide_end_pos = chrom_size - len(prj.config['crop-seq']['cas9'])

    # for each read
    for aln in bam_handle.fetch(reference=chrom + "_chrom"):
        # skip reads
        if (
            aln.is_qcfail or  # failed quality (never happens, but for the future)
            aln.is_secondary or
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
    # Inspect
    fig, axis = plt.subplots(2, sharex=True)
    # number of barcode reads per cell
    sns.distplot(np.log2(1 + reads.groupby(["cell"]).apply(len)), ax=axis[0], kde=False)
    axis[0].set_xlabel("Reads (log2(1 + x))")
    # number of unique barcode reads per cell
    sns.distplot(np.log2(1 + reads.groupby(["cell"])['molecule'].apply(set).apply(len)), ax=axis[1], kde=False)
    axis[1].set_xlabel("Molecules (log2(1 + x))")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.svg"), bbox_inches="tight")

    # process
    u = reads.groupby(
        ['cell', 'molecule', 'chrom'])[
        'distance', 'overlap', 'inside', 'mapping_quality', 'strand_agreeement'].max().reset_index()

    # further reduce molecules by solving chromosome conflicts (assign molecule to chromosome with maximum overlap)
    uu = u.ix[u.groupby(['cell', 'molecule']).apply(lambda x: x['overlap'].argmax())]

    # efficiency of reads in gRNA vs whole construct
    inside_fraction = u.groupby(["cell"])['inside'].sum() / u.groupby(["cell"]).apply(len)
    fig, axis = plt.subplots(1)
    sns.distplot(inside_fraction, bins=20, kde=False)
    axis.set_xlabel("Ratio molecules overlap gRNA / total")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.grna_reads_vs_whole_construct.svg"), bbox_inches="tight")

    # efficiency of reads in gRNA vs whole construct vs number of captured gRNA molecules
    inside_fraction = u.groupby(["cell"])['inside'].sum() / u.groupby(["cell"]).apply(len)

    sns.jointplot(u.groupby(["cell"]).apply(len), u.groupby(["cell"])['inside'].sum())
    axis.set_xlabel("Total molecules in construct per cell")
    axis.set_xlabel("Molecules overlapping gRNA per cell")
    plt.savefig(os.path.join(output_dir, "barcodes_per_cell.all_reads_vs_total_reads.svg"), bbox_inches="tight")
    plt.close("all")

    # remove no overlaps and reads in wrong strand
    u = uu[(uu['overlap'] > 0) & (uu['strand_agreeement'] == 1)]

    # number of unique barcode reads saying inside per cell
    fig, axis = plt.subplots(1)
    sns.distplot(np.log2(1 + u.groupby(["cell"])['molecule'].apply(len)), ax=axis, kde=False)
    axis.set_xlabel("Molecules (log2(1 + x))")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.inside.svg"), bbox_inches="tight")

    # concordance between reads in same cell
    concordant_fraction = 1. / u.groupby(["cell"])['chrom'].nunique()
    fig, axis = plt.subplots(1)
    sns.distplot(concordant_fraction, kde=False)
    axis.set_xlabel("Ratio molecules overlap gRNA / total")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.concordance.svg"), bbox_inches="tight")

    if reads['chrom'].str.contains("Filler_1").any():

        # distribution of reads regarding constructs (for each guide)
        g = sns.FacetGrid(u, col="chrom", sharex=False, sharey=False)
        g.map(sns.distplot, 'distance', kde=False)
        sns.despine(fig)
        g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.distance.svg"), bbox_inches="tight")

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
                    sns.despine(fig)
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
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, "barcodes_per_cell.overlap.in_out.svg"), bbox_inches="tight")


def make_assignment(reads, guide_annotation):
    # Assign
    # unique reads per cell
    # reduce molecules
    u = reads.groupby(
        ['cell', 'molecule', 'chrom'])[
        'distance', 'overlap', 'inside', 'mapping_quality', 'strand_agreeement'].max().reset_index()

    # further reduce molecules by solving chromosome conflicts (assign molecule to chromosome with maximum overlap)
    uu = u.ix[u.groupby(['cell', 'molecule']).apply(lambda x: x['overlap'].argmax())]

    # remove marginal overlaps and reads in wrong strand
    u = uu[(uu['overlap'] > 0) & (uu['strand_agreeement'] == 1)]

    # Get a score (sum of bp covered)
    scores = u.groupby(["cell", 'chrom'])['overlap'].sum()
    scores = scores.reset_index().pivot_table("overlap", "cell", "chrom").fillna(0)

    # assign (get max)
    scores["assignment"] = scores.apply(np.argmax, axis=1)
    scores["score"] = scores[list(set(reads['chrom']))].apply(max, axis=1)
    # give nan to cells with no overlap (this is because argmax pickus a draw)
    scores.loc[scores['score'] == 0, 'assignment'] = pd.np.nan
    scores.loc[scores['score'] == 0, 'score'] = pd.np.nan

    # Get only assigned cells
    scores = scores.dropna()

    # concordance between reads in same cell
    scores['concordance_ratio'] = scores.drop(["assignment", "score"], axis=1).apply(
        lambda x: x / sum(x), axis=1).max(axis=1)

    # Get assigned cells
    assignment = scores.reset_index()[["cell", "assignment", "score", "concordance_ratio"]]

    # Convert to coverage in X times (divide by length of gRNA)
    coverage = scores.drop(['assignment', 'score'], axis=1).apply(
        lambda x: x / float(len(guide_annotation[guide_annotation["oligo_name"] == x.name]["sequence"].squeeze())),
        axis=0)
    coverage["maxscore"] = coverage.apply(max, axis=1)
    coverage["assignment"] = coverage.drop("maxscore", axis=1).apply(np.argmax, axis=1)

    return scores, assignment, coverage


def plot_assignments(scores, assignment, coverage):

    # If number of gRNAs in libarary is less than 20, plot each in a panel separately, else plot all together
    if scores.shape[1] < 22:
        extras = {"col": "assignment", "col_wrap": 4}
    else:
        extras = {}

    # Plot bp covered per cell
    g = sns.FacetGrid(scores[["assignment", "score"]], sharex=False, sharey=False, **extras)
    g.map(sns.distplot, 'score', kde=False)
    g.set(xlim=(0, 1000))
    sns.despine(g.fig)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.bp_covered.svg"), bbox_inches="tight")

    # transform covereage to log
    cov = coverage[["assignment", "maxscore"]]
    cov.loc[:, 'maxscore'] = np.log2(cov['maxscore'])

    try:
        g = sns.FacetGrid(cov, sharex=False, sharey=False, **extras)
        g.map(sns.distplot, 'maxscore', bins=100, kde=False)
        for a in g.axes.flatten():
            a.set_xlabel("log2 coverage")
            sns.despine(g.fig)
        g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.coverage.svg"), bbox_inches="tight")
    except:
        pass

    g = sns.FacetGrid(scores, sharex=True, sharey=False, **extras)
    g.map(sns.distplot, 'concordance_ratio', kde=False)
    sns.despine(g.fig)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.score_concordance_ratio.svg"), bbox_inches="tight")

    # plot assignment stats
    c = assignment['assignment'].value_counts()  # .replace(pd.np.nan, "None")
    fig, axis = plt.subplots(1)
    sns.barplot(c.index, c.values, ax=axis)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=90)
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "barcodes_per_cell.identified.svg"), bbox_inches="tight")

    melted_scores = pd.melt(scores.reset_index(), ['cell', 'assignment', 'score', 'concordance_ratio'])

    g = sns.FacetGrid(melted_scores, sharex=False, sharey=False, **extras)
    g.map(sns.distplot, 'value', bins=200, kde=False)
    sns.despine(fig)
    g.fig.savefig(os.path.join(output_dir, "barcodes_per_cell.scores.distibution.svg"), bbox_inches="tight")

    # calculate abs amount of basepairs overlaping the gRNA of the assigned cell vs all others
    overlap_per_cell = reads.groupby(["cell"])['overlap'].sum()
    overlap_per_guide = reads.groupby(["cell", "chrom"])['overlap'].sum()
    overlap_assignment = scores.apply(lambda x: overlap_per_guide.ix[x.name, x['assignment']] if not pd.isnull(x['assignment']) else pd.np.nan, axis=1)
    overlap_others = overlap_per_cell - overlap_assignment

    sns.jointplot(overlap_others, overlap_assignment, alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.svg"), bbox_inches="tight")
    plt.close("all")
    sns.jointplot(overlap_others, np.log2(1 + overlap_assignment), alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.ylog.svg"), bbox_inches="tight")
    plt.close("all")
    sns.jointplot(np.log2(1 + overlap_others), np.log2(1 + overlap_assignment), alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.bothlog.svg"), bbox_inches="tight")
    plt.close("all")
    sns.jointplot(overlap_others, overlap_assignment, xlim=(-100, overlap_assignment.max() + 100), ylim=(-100, overlap_assignment.max() + 100), alpha=0.1)
    plt.savefig(os.path.join(output_dir, "duplets_assignment_overlap.lims.svg"), bbox_inches="tight")
    plt.close("all")


# Start project, add samples
prj = Project(os.path.join("metadata", "config.yaml"))
prj.add_sample_sheet()
prj.paths.results_dir = results_dir = os.path.join("results")

# get guide annotation
guide_annotation = pd.read_csv(os.path.join("metadata", "guide_annotation.csv"))

for sample in [s for s in prj.samples if hasattr(s, "replicate")]:  # [s for s in prj.samples if hasattr(s, "replicate")]
    output_dir = os.path.join(sample.paths.sample_root, "gRNA_assignment")

    # select gRNAs in respective sample library
    sel_guide_annotation = guide_annotation[guide_annotation['library'] == sample.grna_library]

    # read in alignments
    bam = os.path.join(sample.paths.sample_root, "star_gene_exon_tagged.clean.bam")

    reads = get_reads_in_construct(bam, sel_guide_annotation)
    reads.to_csv(os.path.join(output_dir, "guide_cell_quantification.csv"), index=False)
    reads = pd.read_csv(os.path.join(output_dir, "guide_cell_quantification.csv"))

    # reads in cas9 construct
    cas9_reads = get_reads_in_Cas9_construct(bam)
    cas9_reads.to_csv(os.path.join(output_dir, "cas9_quantification.reads.csv"), index=False)

    cas9_expression = cas9_reads.groupby(['cell'])['molecule'].apply(np.unique).apply(len)
    cas9_expression.reset_index().to_csv(os.path.join(output_dir, "cas9_quantification.counts.csv"), index=False)

    # assign
    scores, assignment, coverage = make_assignment(reads, sel_guide_annotation)
    scores.to_csv(os.path.join(output_dir, "guide_cell_scores.csv"), index=True)
    assignment.to_csv(os.path.join(output_dir, "guide_cell_assignment.csv"), index=False)
    coverage.to_csv(os.path.join(output_dir, "guide_cell_coverage.csv"), index=True)
    scores = pd.read_csv(os.path.join(output_dir, "guide_cell_scores.csv"), index_col=0)
    assignment = pd.read_csv(os.path.join(output_dir, "guide_cell_assignment.csv"))
    coverage = pd.read_csv(os.path.join(output_dir, "guide_cell_coverage.csv"), index_col=0)

    # Plots
    # reads along constructs
    plot_reads_in_constructs(reads)

    # assignment quality/coverage
    plot_assignments(scores, assignment, coverage)

#

#


# Figure 1g
from itertools import chain
colors = sns.color_palette("colorblind")

u6 = prj.config['crop-seq']['u6']
rest = prj.config['crop-seq']['rest']

for sample in [s for s in prj.samples if s.name == "CROP-seq_HEK293T_1_resequenced"]:
    # read in read/construct overlap information
    reads = pd.read_csv(os.path.join(sample.paths.sample_root, "gRNA_assignment", "guide_cell_quantification.csv"))

    # process
    u = reads.groupby(
        ['cell', 'molecule', 'chrom'])[
        'read_start', 'read_end', 'distance', 'overlap', 'inside', 'mapping_quality', 'strand_agreeement'].max().reset_index()
    # further reduce molecules by solving chromosome conflicts (assign molecule to chromosome with maximum overlap)
    uu = u.ix[u.groupby(['cell', 'molecule']).apply(lambda x: x['overlap'].argmax())]
    # remove no overlaps and reads in wrong strand
    u = uu[uu['strand_agreeement'] == 1]

    # select gRNAs in respective sample library
    sel_guide_annotation = guide_annotation[guide_annotation['library'] == sample.grna_library]

    reads2 = u.copy()
    # normalize filler length to match start/end of gRNA
    filler_length = len(sel_guide_annotation.loc[sel_guide_annotation['oligo_name'] == 'Filler_1', 'sequence'].squeeze())
    reads2.loc[
        (reads2["chrom"] == "Filler_1") & (reads2["read_start"] > len(u6) + 20), "read_start"] -= filler_length
    reads2.loc[
        (reads2["chrom"] == "Filler_1") & (reads2["read_end"] > len(u6) + 20), "read_end"] -= filler_length

    # Stacked frequencies of read sequences
    read_data = list()
    for chrom in reads2['chrom'].drop_duplicates():
        if chrom == "Filler_1":
            continue
        read_data.append(
            list(list(chain.from_iterable(
                reads2[(reads2["inside"] == 1) & (reads2["chrom"] == chrom)].apply(
                    lambda x: range(int(x["read_start"]), int(x["read_end"])), axis=1).values))))
    read_data.append(list(list(chain.from_iterable(reads2[reads2["inside"] != 1].apply(lambda x: range(int(x["read_start"]), int(x["read_end"])), axis=1).values))))

    fig, axis = plt.subplots(1, 1, sharex=True)
    axis.hist(
        read_data,
        bins=range(0, len(u6) + 20 + len(rest), 10),
        histtype='barstacked',
        normed=False,
        color=colors[:3] + ["grey"])

    for coor, name in [(0, "startU6"), (len(u6), "start gRNA"), (len(u6) + 20, "start backbone"), (len(u6) + 20 + len(rest), "start polyA")]:
        axis.axvline(coor, 0, 1, linewidth=3, color="black", linestyle="--")
        axis.text(coor, 0.01, name)
    axis.set_xlim((0, len(u6) + 20 + len(rest) + 50))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "figures", "fig1g.reads.stacked.svg"), bbox_inches="tight")
